classdef LSestimatesVAR < handle
    %LSestimatesVAR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        handleDataSample;
        nLags;
        %% Estimates of of the model
        
        A0hat_n_x_1;   % A0 is the estimated constant
        ALhat_n_x_np;  % AL is the estimated lag polynomial
        SigmaHat;      % Sigma is the estimated covariance matrix
        etahat;        % estimated residuals/forecast errors
        
        %% Estimates of IRF
        C;           % C is the estimated reduced form VMA representation, a long matrix (n  x n(MaxHorizons+1) )
        Ccum;        % Ccum is the estimated cumulative VMA representation, a long matrix (n  x n(MaxHorizons+1) )
        G;          %  G matrix: derivative of vec(C) wrt vec(AL)
        Gcum;       %  Derivatives of Ccum wrt vec(AL).
        
        Omega;      % asymptotic covariance matrix for reduced form coefficients
        Omegainv;   % Block inverse of Omega (assumes homoskedasticity)
        
    end
    
    methods
        function obj = LSestimatesVAR( handleDataSample, nLags)
            obj.handleDataSample = handleDataSample;
            obj.nLags = nLags;
            obj = computeLSestimates(obj);
            obj = computeCovariance(obj); 
            obj = computeVMAandDerivatives(obj);
        end
        function thetaHat = getThetaHat(obj)
            thetaHat =  VecAR.thetaFromALSigma( obj.ALhat_n_x_np, obj.SigmaHat );
        end  
        function ALhat_n_x_np = getAL(obj)
            ALhat_n_x_np = obj.ALhat_n_x_np;
        end
        function Omega = getOmega(obj)
            Omega = obj.Omega;
        end          
    end
    methods (Access = private)
        function obj = computeLSestimates(obj)
            [Y,X] = obj.handleDataSample.getYX(obj.nLags);
            
            slopeEstimates = (Y'*X) * ((X'*X)^(-1));
            obj.A0hat_n_x_1  = slopeEstimates(:,1);
            obj.ALhat_n_x_np = slopeEstimates(:,2:end);
            
            obj.etahat = Y' - slopeEstimates * (X');
            reducedT = size(obj.etahat,2);
            obj.SigmaHat = (obj.etahat * obj.etahat') / reducedT; % Covariance matrix
        end
        function obj = computeCovariance(obj)
            [~,X] = obj.handleDataSample.getYX(obj.nLags);
            [obj.Omega, obj.Omegainv] = computeCovarianceOfTheta( X, obj.etahat );
        end
        function obj = computeVMAandDerivatives(obj)
            [obj.C,obj.Ccum] = VecAR.getVMAfromAL(  obj.ALhat_n_x_np,  configFile.MaxHorizons);
            [obj.G,obj.Gcum] = VecAR.getVMAderivatives(  obj.ALhat_n_x_np,  configFile.MaxHorizons);
        end
    end
    
end 

function [OmegaHat,OmegaHatinv] = computeCovarianceOfTheta(X,eta)
% -------------------------------------------------------------------------
% Computes the asymptotic variance of [vec(Ahat)',vech(Sigmahat)']'
% Please, refer to  Lütkepohl H. New introduction to multiple time series analysis. ? Springer, 2007.
%
% Inputs:
% - X: matrix of dimension (T times n) containing the time series
% - eta: reduced-form residuals
% Outputs:
% - OmegaHat: asymptotic variance of [vec(Ahat)',vech(Sigmahat)']'
%
% This version: March 23, 2017
% Last edited by Bulat Gafarov
% -------------------------------------------------------------------------

n = size(eta,1);
n_x_p = size(X,2) - 1;
p = n_x_p / n;

scedasticity = configFile.scedasticity;
vecFromVech = converter.getVecFromVech(n);

%% suppress default warning for near singularity
warning('off','MATLAB:nearlySingularMatrix')



switch scedasticity
    case 'homo'
        
        %% Definitions
        XSVARp = X(:,2:end);
        T1aux = size(eta,2); %This is the number of time periods
        Sigmau = (eta*eta')/(size(eta,2));
        Gamma = XSVARp'*XSVARp/T1aux;
        DK = vecFromVech;     % duplication matrix such that vec(Sigma)=D*vech(Sigma)
        DKplus = (DK'*DK)\DK';
        
        
        %% Construct Sigma_a (p.118 in Luetkepohl, 1990)
        Sigmaa = kron(Gamma\eye(size(Gamma)),Sigmau);
        na = size(Sigmaa,1);
        
        
        %% Construct Sigma_s (p.118 in Luetkepohl, 1990)
        Sigmas = 2*DKplus*kron(Sigmau,Sigmau)*DKplus';
        ns = size(Sigmas,1);
        
        
        %% Omega
        OmegaHat = [Sigmaa, zeros(na,ns); zeros(ns,na), Sigmas];
        OmegaHatinv = OmegaHat\eye(size(OmegaHat));
        
        
    case 'hetero'
        
        %% Definitions
         
        XSVARp = X;
        matagg = [XSVARp,eta']'; %The columns of this vector are (1;X_t; eta_t)
        T1aux = size(eta,2);    %This is the number of time periods
        T2aux = size(matagg,1); %This is the column dimension of (1;X_t;eta_t)
        
        
        %%
        etaaux = reshape(eta,[n,1,T1aux]); %Each 2-D page contains eta_t
        mataggaux = permute(reshape(matagg,[T2aux,1,T1aux]),[2,1,3]); %Each 2-D page contains (1,X_t',\eta_t')
        auxeta = bsxfun(@plus,bsxfun(@times,etaaux,mataggaux),-mean(bsxfun(@times,etaaux,mataggaux),3));
        %Each 2-D page contains [eta_t, eta_t X_t', eta_t eta_t'-Sigma];
        vecAss1= reshape(auxeta,[n+(p*(n^2))+n^2,1,T1aux]);
        %Each 2-D page contains [eta_t; vec(eta_tX_t') ; vec(eta_t*eta_t'-Sigma) ]
        WhatAss1 = sum(bsxfun(@times,vecAss1,permute(vecAss1,[2,1,3])),3)./T1aux;
        %This is the covariance matrix we are interested in
        
        
        %% Construct the selector matrix VauxInverse that gives: vech(Sigma)=Vaux*vec(Sigma)
        I = eye(n);
        VauxInverse = kron(I(1,:),I);
        for i=2:n
            VauxInverse = [VauxInverse; kron(I(i,:),I(i:end,:))];
        end
        
        
        %% Check bad conditioning
        test = n^2*p + n*(n+1)/2 < T1aux;
        if false(test)
            error('Warning: The covariance matrix is not invertible because n^2p+n(n+1)/2<T!')
        end
        
        
        %% This is the estimator for matrix W in Assumption 1
        Mhat = [kron([zeros(n*p,1),eye(n*p)]*(XSVARp'*XSVARp./T1aux)^(-1),eye(n)), zeros((n^2)*p,n^2); ...
            zeros(n*(n+1)/2,((n^2)*p)+n), VauxInverse];
        OmegaHat = (Mhat)*(WhatAss1)*(Mhat');
        OmegaHatinv = OmegaHat\eye(size(OmegaHat));
    otherwise
        error('Warning: choose homo or heteorscedasticity option')
end


%% Warning for near-singularity
test = rcond(OmegaHat);
if test<eps
    fprintf('... warning: OmegaHat covariance matrix is close to singular!\n')
end
warning('on')
end

