classdef (Abstract) VecAR
    %VECAR Summary of this class goes here
    %   Detailed explanation goes here
    properties  ( Access = protected)
        nLags = [] ;  %% The number "p" of lags in the SVAR model
        scedasticity;
        config;     % handle to an object of SVARconfiguration class
        cache;
    end
 
    methods
        function IRFObjectiveFunctions = getIRFObjectiveFunctions(obj)
            switch (obj.config.isCumulativeIRF)
                case 'yes'
                    reducedFormIRFnonCum = obj.getVMA_ts_sh_ho;
                    IRFObjectiveFunctions = cumsum(reducedFormIRFnonCum,3);
                otherwise
                    IRFObjectiveFunctions = obj.getVMA_ts_sh_ho;
            end;
            
        end
        function IRFObjectiveFunctions = getIRFObjectiveFunctionsDerivatives(obj)
            switch (obj.config.isCumulativeIRF)
                case 'yes'
                    reducedFormIRFnonCum = obj.getVMADerivatives_ts_sh_ho_dAL;
                    IRFObjectiveFunctions = cumsum(reducedFormIRFnonCum,3);
                otherwise
                    IRFObjectiveFunctions = obj.getVMADerivatives_ts_sh_ho_dAL;
            end;
            
        end
        function configHandle = getConfig (obj)
            configHandle = obj.config;
        end
        function obj = precomputeCache(obj)
            obj.cache.vecFromVech = converter.getVecFromVech(obj.getN);
            obj.cache.G = getVMADerivatives_ts_sh_ho_dAL(obj)  ;
            obj.cache.VMA_ts_sh_ho = getVMA_ts_sh_ho(obj);
            obj.cache.Sigma = getSigma(obj);
        end
        function [AL_n_x_np,Sigma] =  ALSigmaFromThetaNandP(obj, theta, nLags) 
            nShocks = obj.getN;
            dAL = nShocks*nShocks*nLags;
            vecAL = theta(1:dAL);
            AL_n_x_np  = reshape(vecAL,[nShocks, nShocks*nLags] );
            
            vechSigma = theta((dAL+1):end);
            
            if isfield( obj.cache, 'vecFromVech')
                vecFromVech = obj.cache.vecFromVech;
            else
                vecFromVech = converter.getVecFromVech( nShocks);
            end
            vecSigma =  vecFromVech * vechSigma;
            Sigma = reshape(vecSigma,[nShocks,nShocks]);
        end
        function MaxHorizons = getMaxHorizons(obj)
            MaxHorizons = obj.config.nNoncontemoraneousHorizons+1;
        end
    end
    methods (Abstract)
        getVMA_ts_sh_ho(obj);
        getVMADerivatives_ts_sh_ho_dAL(obj);
        getSigma(obj);
        getN(obj);
        getTheta(obj);
        getNames(obj);
        getUnitsOfMeasurement(obj);
        getCovarianceOfThetaT(obj);
    end
   
    methods (Access = public, Static)
        function theta = thetaFromALSigma(AL_n_x_np,Sigma)
            %% this funciton computes the vector of the reduced form parameters, theta vector
            [n,p] = VecAR.getNPfromAL(AL_n_x_np);
            
            vecAL = reshape(AL_n_x_np,  [ (n^2)*p ,1]);
            vecSigma = reshape( Sigma, [n^2,1]);
            
            vechSigmaFromVec = converter.getVechFromVec(n);
            vechSigma = vechSigmaFromVec * vecSigma;
            
            theta = [vecAL; vechSigma];
        end
        function VMA_ts_sh_ho  = getVMAfromAL( AL_n_x_np, nNoncontemoraneousHorizons)
            % -------------------------------------------------------------------------
            % Transforms the A(L) parameters of a reduced-form VAR
            % into the coefficients C of the MA representation.
            %
            % Inputs:
            % - AL: VAR model coefficients
            % - nNoncontemoraneousHorizons : forecast horizon
            % Outputs:
            % - C: MA representation coefficients
            %
            % This version: March 31, 2015
            % -------------------------------------------------------------------------
            
            
            %% Reshape AL into a 3-D array
            [n,p] = VecAR.getNPfromAL(AL_n_x_np);
            AL_n_x_n_x_p = reshape(AL_n_x_np,[n,n,p]);
            
            %% Initialize the value of the auxiliary array vecALrevT
            vecALrevT = zeros(n,n,nNoncontemoraneousHorizons);
            for i=1:nNoncontemoraneousHorizons
                if i<(nNoncontemoraneousHorizons-p)+1
                    vecALrevT(:,:,i) = zeros(n,n);
                else
                    vecALrevT(:,:,i) = AL_n_x_n_x_p(:,:,(nNoncontemoraneousHorizons-i)+1)';
                end
            end
            vecALrevT = reshape(vecALrevT,[n,n*nNoncontemoraneousHorizons]);
            
            
            %% MA coefficients
            C = repmat(AL_n_x_n_x_p(:,:,1),[1,nNoncontemoraneousHorizons]);
            for i=1:nNoncontemoraneousHorizons-1
                C(:,(n*i)+1:(n*(i+1))) = [eye(n),C(:,1:n*i)] * vecALrevT(:,(nNoncontemoraneousHorizons * n-(n*(i+1)))+1:end)';
            end

            
            %% cumulative MA coefficients
            VMA = [eye(n), C];
            VMA_ts_sh_ho = reshape(VMA,[n,n,(nNoncontemoraneousHorizons+1)]);
        end
        function G_ts_sh_ho_dAL = getVMAderivatives(AL_n_x_np, hori)
            % -------------------------------------------------------------------------
            % Computes the derivatives of vec(C) wrt vec(A) based on
            % Lütkepohl H. New introduction to multiple time series analysis. ? Springer, 2007.
            %
            % Inputs:
            % - AL: VAR model coefficients
            % - hori: forecast horizon
            % Outputs:
            % - G: derivatives of C wrt A
            %
            % This version: March 23, 2015
            % -------------------------------------------------------------------------
            
            VMA_ts_sh_ho = VecAR.getVMAfromAL( AL_n_x_np, hori);
            [n,p] = VecAR.getNPfromAL( AL_n_x_np );

            %% A and J matrices in Lutkepohl's formula for the derivative of C with respect to A
            J = [eye(n), zeros(n,(p-1)*n)];
            Alutkepohl = [AL_n_x_np; eye(n*(p-1)), zeros(n*(p-1),n)];
            
            %% AJ is a 3D array that contains A^(k-1) J' in the kth 2D page of the the 3D array
            AJ = zeros(n*p, n, hori);
            for k=1:hori
                AJ(:,:,k) = ((Alutkepohl)^(k-1)) * J';
            end
            
            %% matrix [ JA'^0; JA'^1; ... J'A^{k-1} ];
            JAp = reshape(AJ, [n*p,n*hori])';
            %% G matrices
            AJaux = zeros(size(JAp,1)*n, size(JAp,2)*n, hori);
            for i=1:hori
                AJaux(((n^2)*(i-1))+1:end,:,i) = kron(JAp(1:n*(hori+1-i),:), VMA_ts_sh_ho(:,:,i));
            end
            Gaux = permute(reshape(sum(AJaux,3)', [(n^2)*p, n^2, hori]), [2,1,3]);
            G_tsxsh_dAL_ho = zeros(size(Gaux,1), size(Gaux,2), size(Gaux,3)+1);
            G_tsxsh_dAL_ho(:,:,2:end) = Gaux;
            
            G_dAL_ho_tsxsh = permute(G_tsxsh_dAL_ho,[2,3,1]);
            dAL = size(Gaux,2);
            MaxHorizon = hori+1;
            G_dAL_ho_ts_sh = reshape(G_dAL_ho_tsxsh,[dAL,MaxHorizon,n,n]);
            G_ts_sh_ho_dAL = permute(G_dAL_ho_ts_sh,[3,4,2,1] );
        end
        function [n,p] = getNPfromAL(AL_n_x_np)
            n = size(AL_n_x_np,1);
            n_x_p = size(AL_n_x_np, 2);
            p = n_x_p / n;
        end
        function stationarity(AL)
            % -------------------------------------------------------------------------
            % Checks whether reduced-form VAR model is covariance stationary
            %
            % Inputs:
            % - AL: VAR model coefficients
            % Outputs: (none)
            %
            % This version: March 31, 2015
            % -------------------------------------------------------------------------
            [n,p] = VecAR.getNPfromAL(AL);
            A = [ AL; eye( n*(p-1) ), zeros( n*(p-1), n)];
            
            eigenvaluesOfA = eigs(A, size(A,1)-2); % FIXME: WHY -2?
            
            if max(abs(eigenvaluesOfA))>= 1;
                disp('Warning: VAR model is not covariance stationary!');
            else
                disp('VAR model is covariance stationary');
            end
            
            
        end

        
    end
    
end

