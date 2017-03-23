classdef LSestimatesVAR < handle
    %LSestimatesVAR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private) 
        A0hat_n_x_1;   % A0 is the estimated constant
        ALhat_n_x_np;  % AL is the estimated lag polynomial
        SigmaHat;      % Sigma is the estimated covariance matrix
        etahat;        % estimated residuals/forecast errors
        handleDataSample;
        nLags;
    end
    
    methods
        function obj = LSestimatesVAR( handleDataSample, nLags)
            obj.handleDataSample = handleDataSample;
            obj.nLags = nLags;
            
            [Y,X] = handleDataSample.getYX(obj.nLags);
            
            slopeEstimates = (Y'*X) * ((X'*X)^(-1));
            obj.A0hat_n_x_1  = slopeEstimates(:,1);
            obj.ALhat_n_x_np = slopeEstimates(:,2:end);
            
            obj.etahat = Y' - slopeEstimates * (X');
            reducedT = size(obj.etahat,2);
            obj.SigmaHat = (obj.etahat * obj.etahat') / reducedT; % Covariance matrix
        end
        function thetaHat = getThetaHat(obj)
            thetaHat =  ALSigmaToTheta( obj.ALhat_n_x_np, obj.SigmaHat );
        end
        function obj = computeCovariance(obj)
            %% Asymptotic covariance matrix for reduced form coefficients
            [obj.Omega,obj.Omegainv] = computeCovarianceOfTheta(obj.nLags,obj.X,obj.etahat,obj.vecFromVech,obj.config.scedasticity);
        end
    end
    
end



function theta = ALSigmaToTheta(AL_n_x_np,Sigma)
%% this funciton computes the vector of the reduced form parameters, theta vector
n = size(Sigma,1);
n_x_p = size(AL_n_x_np, 2);
p = n_x_p / n;

vecAL = reshape(AL_n_x_np,[(n^2)*p,1]);
vecSigma = reshape( Sigma, [n^2,1]);

vechSigmaFromVec = converter.getVechFromVec(n);
vechSigma = vechSigmaFromVec * vecSigma;

theta = [vecAL; vechSigma];

end

 
function [OmegaHat,OmegaHatinv] = computeCovarianceOfTheta(p,X,eta,vecFromVech,scedasticity)
% -------------------------------------------------------------------------
% Computes the asymptotic variance of [vec(Ahat)',vech(Sigmahat)']'
% Please, refer to  Lütkepohl H. New introduction to multiple time series analysis. ? Springer, 2007.
%
% Inputs:
% - p: lag order
% - Sigma: covariance matrix of residuals
% - Y: matrix of dimension (T times n) containing the time series
% - eta: reduced-form residuals
% Outputs:
% - OmegaHat: asymptotic variance of [vec(Ahat)',vech(Sigmahat)']'
%
% This version: March 21, 2017
% Last edited by Bulat Gafarov
% -------------------------------------------------------------------------


%% suppress default warning for near singularity
warning('off','MATLAB:nearlySingularMatrix')



switch scedasticity
    case 'homo'
        
        %% Definitions
        n = size(eta,1);
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
        n = size(eta,1);
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
disp('  Asym. covariance for (AL,Sigma) is computed successfully.')
end