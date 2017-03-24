classdef VecAR
    %VECAR Summary of this class goes here
    %   Detailed explanation goes here
    properties (Access = protected)
        nLags = [] ;  %% The number "p" of lags in the SVAR model
        config;     % handle to an object of configFile class
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
        function [C,Ccum] = getVMAfromAL( AL_n_x_np, hori)
            % -------------------------------------------------------------------------
            % Transforms the A(L) parameters of a reduced-form VAR
            % into the coefficients C of the MA representation.
            %
            % Inputs:
            % - AL: VAR model coefficients
            % - hori: forecast horizon
            % Outputs:
            % - C: MA representation coefficients
            %
            % This version: March 31, 2015
            % -------------------------------------------------------------------------
            
            
            %% Reshape AL into a 3-D array
            [n,p] = VecAR.getNPfromAL(AL_n_x_np);
            AL_n_x_n_x_p = reshape(AL_n_x_np,[n,n,p]);
            
            %% Initialize the value of the auxiliary array vecALrevT
            vecALrevT = zeros(n,n,hori);
            for i=1:hori
                if i<(hori-p)+1
                    vecALrevT(:,:,i) = zeros(n,n);
                else
                    vecALrevT(:,:,i) = AL_n_x_n_x_p(:,:,(hori-i)+1)';
                end
            end
            vecALrevT = reshape(vecALrevT,[n,n*hori]);
            
            
            %% MA coefficients
            C = repmat(AL_n_x_n_x_p(:,:,1),[1,hori]);
            for i=1:hori-1
                C(:,(n*i)+1:(n*(i+1))) = [eye(n),C(:,1:n*i)] * vecALrevT(:,(hori*n-(n*(i+1)))+1:end)';
            end
            
            
            %% cumulative MA coefficients
            Ctmp = [eye(n), C];
            Chataux = cumsum(reshape(Ctmp,[n,n,(hori+1)]),3);
            Ccum = reshape(Chataux,[n,n*(hori+1)]);
            Ccum = Ccum(:,(n+1):end);
            
        end
        function [G,Gcum] = getVMAderivatives(AL_n_x_np, hori)
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
            
            [C,~] = VecAR.getVMAfromAL( AL_n_x_np, hori);
            [n,p] = VecAR.getNPfromAL( AL_n_x_np );

            %% A and J matrices in Lutkepohl's formula for the derivative of C with respect to A
            J = [eye(n), zeros(n,(p-1)*n)];
            Alutkepohl = [AL_n_x_np; eye(n*(p-1)),zeros(n*(p-1),n)];
            
            
            %% AJ is a 3D array that contains A^(k-1) J' in the kth 2D page of the the 3D array
            AJ = zeros(n*p, n, hori);
            for k=1:hori
                AJ(:,:,k) = ((Alutkepohl)^(k-1)) * J';
            end
            
            %% matrix [ JA'^0; JA'^1; ... J'A^{k-1} ];
            JAp = reshape(AJ, [n*p,n*hori])';
            
            
            %% G matrices
            AJaux = zeros(size(JAp,1)*n, size(JAp,2)*n, hori);
            Caux = reshape([eye(n), C(:,1:(hori-1)*n)], [n,n,hori]);
            for i=1:hori
                AJaux(((n^2)*(i-1))+1:end,:,i) = kron(JAp(1:n*(hori+1-i),:), Caux(:,:,i));
            end
            Gaux = permute(reshape(sum(AJaux,3)', [(n^2)*p, n^2, hori]), [2,1,3]);
            G = zeros(size(Gaux,1), size(Gaux,2), size(Gaux,3)+1);
            G(:,:,2:end) = Gaux;
            
            
            %% Cumulative version of G matrices
            Gcum = cumsum(G,3);
            
            
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
        
        
        
        function simVAR   = simulateVAR(AL_n_x_np,T,nMC,etahat)
            % -------------------------------------------------------------------------
            % This function simulate nMC samples of length T based on VAR model with lag polynomial AL
            % and covariance matrix Sigma
            %
            % Inputs:
            % - AL: VAR model coefficients
            % - T: timer series length
            % - nMC: number of bootstrap samples
            % - eta: VAR model residuals
            % Outputs:
            % - simVAR: (T x n x nMC) 3d array of simulated data
            %
            % This version: February 24, 2015
            % -------------------------------------------------------------------------
            
            % To do list
            % - use MA representation to vectorize the simulation for better speed
            
            [n,np] = size(AL_n_x_np);
            lags = np/n;
            burnIn = 1000;
            etahat =  etahat -mean(etahat,2)* ones( 1, T-lags) ;
            
            eta = etahat( randi([1,(T-lags)],nMC*T+burnIn,n) );
            TSLtemp =  zeros( nMC * T + burnIn, n);
            for iT = ( lags+1):( nMC * T + burnIn)
                TSLtemp(iT,:) =  ( AL_n_x_np * reshape( (TSLtemp((iT-1):-1:(iT- lags),:))',[ n *  lags,1]) )' + eta(iT,:);
            end
            
            simVAR = permute ( reshape(TSLtemp((burnIn+1):( nMC* T+burnIn),:),[ T ,  nMC ,  n ]),   [1 3 2]);
            
        end


    end
    
end

