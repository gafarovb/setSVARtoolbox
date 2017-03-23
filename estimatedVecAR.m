classdef estimatedVecAR 
    %VecAR is a class for a time series collection
    %   Detailed explanation goes here

    %%
    properties (Access = private)
        %% seeds for random numbers
        seedMC = NaN; % seedMC: seed for a given MC sample. Use NaN for original dataset
        %% Data
        dataSample;
    end
    %%
    properties (SetAccess = protected)   
        config; % object of configFile class
         %% model characteristics
        nLags = 0 ;  %% The number "p" of lags in the SVAR model
        T = 0 ;       % Number of time periods
%         d = 0 ;          % Length of vector of reduced-form parameters
        estimates;

%         %% Variance
%         Omega;      % asymptotic covariance matrix for reduced form coefficients
%         Omegainv;   % Block inverse of Omega (assumes homoskedasticity)
        %% Estimates of IRF
        C;           % C is the estimated reduced form VMA representation, a long matrix (n  x n(MaxHorizons+1) )
        Ccum;        % Ccum is the estimated cumulative VMA representation, a long matrix (n  x n(MaxHorizons+1) )
        %% todo : check if these properties are necessary
        %% Auxilary properties
        G;          % Derivatives of C wrt AL.
        Gcum;       %  Derivatives of Ccum wrt AL.
        bootCSigma; % Bootstrap after boostrap samples of of C wrt AL.
    end
    %% ********************* Methods ***************************************
    methods
        function obj = VecAR

            obj.config = configFile;
            obj.nLags = obj.config.nLags ;
            obj.dataSample = multivariateTimeSeries;  
            obj.estimates = LSestimatesVAR( obj.dataSample, obj.nLags);
%             obj = computeCovariance(obj);
%             obj = VMArepresentationAndDerivatives(obj);            
            
        end

        function obj = VMArepresentationAndDerivatives(obj)
            %% VMA coefficients
            [obj.C,obj.Ccum] = computeVMArepresentation(obj.AL,obj.nLags,obj.config.MaxHorizons);
            %% G matrix: derivative of vec(C) wrt vec(AL)
            [obj.G,obj.Gcum] = computeG_VMAderivatives(obj.AL,obj.C,obj.nLags,obj.config.MaxHorizons,obj.getN);
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
            stationarity(obj.AL,obj.nLags,obj.getN)
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




%% functions from folder funcRForm

function [C,Ccum] = computeVMArepresentation(AL,p,hori)
% -------------------------------------------------------------------------
% Transforms the A(L) parameters of a reduced-form VAR
% into the coefficients C of the MA representation.
%
% Inputs:
% - AL: VAR model coefficients
% - p: number of lags in the VAR model
% - hori: forecast horizon
% Outputs:
% - C: MA representation coefficients
%
% This version: March 31, 2015
% -------------------------------------------------------------------------


%% Reshape AL into a 3-D array
n = size(AL,1);
vecAL = reshape(AL,[n,n,p]);

%% Initialize the value of the auxiliary array vecALrevT
vecALrevT = zeros(n,n,hori);
for i=1:hori
    if i<(hori-p)+1
        vecALrevT(:,:,i) = zeros(n,n);
    else
        vecALrevT(:,:,i) = vecAL(:,:,(hori-i)+1)';
    end
end
vecALrevT = reshape(vecALrevT,[n,n*hori]);


%% MA coefficients
C = repmat(vecAL(:,:,1),[1,hori]);
for i=1:hori-1
    C(:,(n*i)+1:(n*(i+1))) = [eye(n),C(:,1:n*i)] * vecALrevT(:,(hori*n-(n*(i+1)))+1:end)';
end


%% cumulative MA coefficients
Ctmp = [eye(n), C];
Chataux = cumsum(reshape(Ctmp,[n,n,(hori+1)]),3);
Ccum = reshape(Chataux,[n,n*(hori+1)]);
Ccum = Ccum(:,(n+1):end);

end
function [G,Gcum] = computeG_VMAderivatives(AL,C,p,hori,n)
% -------------------------------------------------------------------------
% Computes the derivatives of vec(C) wrt vec(A) based on
% Lütkepohl H. New introduction to multiple time series analysis. ? Springer, 2007.
%
% Inputs:
% - AL: VAR model coefficients
% - C: MA representation coefficients
% - p: lag order
% - hori: forecast horizon
% - n: number of variables
% Outputs:
% - G: derivatives of C wrt A
%
% This version: March 31, 2015
% -------------------------------------------------------------------------


%% A and J matrices in Lutkepohl's formula for the derivative of C with respect to A
J = [eye(n), zeros(n,(p-1)*n)];
Alut = [AL; eye(n*(p-1)),zeros(n*(p-1),n)];


%% AJ is a 3D array that contains A^(k-1) J' in the kth 2D page of the the 3D array
AJ = zeros(n*p, n, hori);
for k=1:hori
    AJ(:,:,k) = ((Alut)^(k-1)) * J';
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
function simVAR   = simulateVAR(AL,T,nMC,etahat)
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

[n,np] = size(AL);
lags = np/n;
burnIn = 1000;
etahat =  etahat -mean(etahat,2)* ones( 1, T-lags) ;

eta = etahat( randi([1,(T-lags)],nMC*T+burnIn,n) );
TSLtemp =  zeros( nMC * T + burnIn, n);
for iT = ( lags+1):( nMC * T + burnIn)
    TSLtemp(iT,:) =  ( AL * reshape( (TSLtemp((iT-1):-1:(iT- lags),:))',[ n *  lags,1]) )' + eta(iT,:);
end

simVAR = permute ( reshape(TSLtemp((burnIn+1):( nMC* T+burnIn),:),[ T ,  nMC ,  n ]),   [1 3 2]);

end
function stationarity(AL,p,n)
% -------------------------------------------------------------------------
% Checks whether reduced-form VAR model is covariance stationary
%
% Inputs:
% - AL: VAR model coefficients
% - p: number of lags in the VAR model
% - n: length of time series dimension
% Outputs: (none)
%
% This version: March 31, 2015
% -------------------------------------------------------------------------


%% Definitions
A = [AL; eye(n*(p-1)),zeros(n*(p-1),n)];


%% Eigenvalues
e = eigs(A,size(A,1)-2);
r = real(e);
i = imag(e);
dd = sqrt(r.^2+i.^2);


if max(dd)>=1;
    disp('Warning: VAR model is not covariance stationary!');
else
    disp('VAR model is covariance stationary');
end


end

%********************************************************
%********************************************************
%********************************************************
%********************************************************
%********************************************************