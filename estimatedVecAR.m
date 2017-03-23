classdef estimatedVecAR < VecAR
    %VecAR is a class for a time series collection
    %   Detailed explanation goes here

    %%
    properties (Access = private)
        seedMC = NaN; % seedMC: seed for a given MC sample. Use NaN for original dataset
        dataSample;
    end
    
    properties (SetAccess = protected)   
        config ;     % object of configFile class
        nLags = 0 ;  %% The number "p" of lags in the SVAR model
        estimates;
        %% todo : check if these properties are necessary
        bootCSigma; % Bootstrap after boostrap samples of of C wrt AL.
    end
    
    
    methods
        function obj = estimatedVecAR
            obj.config = configFile;
            obj.nLags = obj.config.nLags ;
            obj.dataSample = multivariateTimeSeries;  
            obj.estimates = LSestimatesVAR( obj.dataSample, obj.nLags);
        end
        function obj = eraseData(obj)
            % this method erases information about the sample but leaves
            % the dimensions of the data intact
            obj.tsInColumns = obj.tsInColumns * NaN;
            obj.etahat = obj.etahat * NaN;
            obj.X =  obj.X * NaN;
            obj.Y =  obj.Y * NaN;
        end
        function stationarityTest(obj)
            % -------------------------------------------------------------------------
            % Checks whether reduced-form VAR model is covariance stationary
            % -------------------------------------------------------------------------
            VecAR.stationarity(obj.estimates.getAL)
        end
        function obj = bab(obj)
            % -------------------------------------------------------------------------
            % bootstrap after bootstrap
            % See Kilian (1998) ReStat
            %
            % Inputs:
            % - obj.tsInColumns  = series: matrix of dimension (T times n) containing the time series
            % - obj.nLags = p: model lag order
            % - obj.config.MaxHorizons = hori: model horizon
            % - obj.config.babn1 = n1: number of bootstrap samples before bias correction
            % - obj.config.babn2 = n2: number of bootstrap samples after  bias correction
            % Outputs:
            % - output: structure containing all reduced-form model objects
            %           after bootstrapping
            %
            % This version: March 31, 2015
            % -------------------------------------------------------------------------
            
            
            %%  read inputs from the VecAR object
            series = obj.tsInColumns;
            p      = obj.nLags;
            hori   = obj.config.MaxHorizons ;
            n1 = obj.config.babn1;
            n2 = obj.config.babn2;
            AL = obj.AL;
            sigma = obj.Sigma;
            etahat = obj.etahat;
            T = obj.T;
            
            %% First stage bootrap
            simVAR = simulateVAR(AL,T,n1,etahat);
            [n,np] = size(AL);
            p = np/n;
            ALar = zeros(n,np,n1);
            sigmaAr = zeros(n,n,n1);
            for iMC = 1:n1
                seriesMC = simVAR(:,:,iMC);
                seriesMC = seriesMC - ones(T,1)* mean(seriesMC);
                [ALar(:,:,iMC),sigmaAr(:,:,iMC),~] = computeLSestimatesOfVAR(seriesMC,p);
                % stationarity checks
                stationarity(ALar(:,:,iMC),p,n)
            end
            
            
            %% Use resampled data to correct bias
            meanALb = mean(ALar,3);
            ALbc = 2*AL-meanALb ;
            sigmabc = 2*sigma - mean(sigmaAr,3);
            companion = [ALbc;eye(n*(p-1),np)] ;
            
            
            %% make it a separate function
            if max(abs(eig(companion)))>0.95
                deltaU = 1;
                deltaL = 0;
                delta  = 0.5;
                while abs(delta-deltaU)+abs(delta-deltaL)>0.01
                    ALbc = AL + delta*(AL-meanALb);
                    companion = [ALbc;eye(n*(p-1),np)] ;
                    if (max(abs(eig(companion)))>0.95)
                        deltaU = delta;
                        delta = (deltaU + deltaL)/2 ;
                    else
                        deltaL = delta;
                        delta = (deltaU + deltaL)/2 ;
                    end
                end
            end
            
            
            %% Get corrected residuals
            aux = lagmatrix(series,[1:1:p]);
            XSVARp = bsxfun(@plus,aux((p+1):end,:),-mean(aux((p+1):end,:),1));
            etahat = series((p+1):end,:)'-ALbc*(XSVARp');
            
            
            %% Resample data
            simVAR = simulateVAR(ALbc,T,n2,etahat);
            CAr = zeros(n ,hori*n , n2 );
            CcumAr = zeros(n ,hori*n , n2 );
            ALar = zeros(n ,np , n2 );
            sigmaAr = zeros(n ,n , n2 );
            for iMC = 1:n2
                seriesMC = simVAR(:,:,iMC);
                seriesMC = seriesMC - ones(T,1)* mean(seriesMC);
                [ALar(:,:,iMC),sigmaAr(:,:,iMC),~] = computeLSestimatesOfVAR(seriesMC,p);
                [CAr(:,:,iMC)]=computeVMArepresentation(ALar(:,:,iMC),p,hori);
                Ctmp = [eye(n),CAr(:,:,iMC)];
                Chataux = cumsum(reshape(Ctmp,[n,n,(hori+1)]),3);
                Ccum = reshape(Chataux,[n,n*(hori+1)]);
                CcumAr(:,:,iMC) = Ccum(:,(n+1):end);
            end
            
            
            %% save output
            obj.bootCSigma = struct(...
                'C',CAr,...
                'Sigma',sigmaAr,...
                'Ccum',CcumAr,...
                'ALbc',ALbc,...
                'Sigmabc',sigmabc);
            
        end
        function [bic,aic,hqic] = bicaic(obj)
            % -------------------------------------------------------------------------
            % This function selects number of lags in the VAR model for "series"
            % based on AIC,BIC and HQ information criteria
            %
            % Inputs:
            % - series: matrix of dimension (T times n) containing the time series
            % Outputs:
            % - optimal lag-length by information criterion
            %
            % This version: March 31, 2015
            % Please, cite Gafarov, B. and Montiel-Olea, J.L. (2015)
            % "ON THE MAXIMUM AND MINIMUM RESPONSE TO AN IMPULSE IN SVARS"
            % -------------------------------------------------------------------------
            series= obj.tsInColumns;
            pmax = 24; % - pmax: maximum possible number of lags
            
            %% Definitions
            [T,N]=size(series(pmax+1:end,:)); % consider first pmax periods as presample
            
            
            %% Allocate memory
            aic = zeros(pmax,1);
            bic = zeros(pmax,1);
            hqic = zeros(pmax,1);
            
            
            %% Compute information criteria
            for px = 1:pmax
                series_ = series(pmax-px+1:end,:);
                %     series_ = series;
                [~,~,eta,~] = computeLSestimatesOfVAR(series_,px); %Apply the function RForm_VAR.m to estimate reduced form parameters
                lhd = log(det((eta*eta')/T));
                pty = px*N^2/T;
                aic(px)  = lhd + 2*pty;
                bic(px)  = lhd + log(T)*pty;
                hqic(px) = lhd + 2*log(log(T))*pty;
            end
            
            
            %% Optimal lag
            [~,aic] = min(aic);
            [~,bic] = min(bic);
            [~,hqic] = min(hqic);
            
        end
        function objSimulated = resampleTheta(obj,seedMC)
            if nargin>1
                %%  simulates theta from the asymptotic distribution
                %   and erases the data from objSimulated;
                
                objSimulated = eraseData(obj);
                
                objSimulated.seedMC = seedMC;
                rng('default');
                rng(objSimulated.seedMC,'twister');
                errors = randn(objSimulated.d,1);
                
                objSimulated.theta  =  obj.theta+((obj.Omega)^(1/2)/(obj.T^.5))*errors;
                objSimulated.AL     = reshape(objSimulated.theta(1:(obj.getN^2)*obj.nLags),[obj.getN,obj.getN*obj.nLags]);
                objSimulated.Sigma  = reshape(obj.vecFromVech*objSimulated.theta((obj.getN^2)*obj.nLags+1:end,1),[obj.getN,obj.getN]);
                
                objSimulated = VMArepresentation(objSimulated);
                
            else
                disp('Error: Provide a seed to resample theta.')
            end
            
        end
        function n = getN(obj)
            n = obj.dataSample.countTS;
        end
        function d = countParameters(obj)
            n = obj.getN ;
            p = obj.nLags;
            d = n * n * p + n * (n+1) / 2;   
        end
    end
    
end




