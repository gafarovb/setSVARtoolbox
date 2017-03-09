classdef VecAR 
    %VAR is a class for a time series collection 
    %   Detailed explanation goes here
 
    properties (Access = public) 
      names = []; %% Labels for TS         
    end
    properties (Constant)      % todo: make it a cofiguration file
        MaxHorizons = 23;      % the number of horizons to compute the IRF 
        hetscedOmega = 0;      % 0=homoscedastic Omega, 1=heteroscedastic Omega
        babn1 = 1000;          % bab: number of bootstrap samples before bias correction
	    babn2 = 1000;          % bab: number of bootstrap samples after  bias correction
    end
  
%%  
    properties (Access = private)
    %% Data  
      X;           % Regressors
      Y;           % Regressand 
    end
%%    
    properties (Access = protected) 
  %%  Data         
      data = [] ;      
  %% model characteristics
      nLags = 0 ;  %% The number "p" of lags in the SVAR model
      n = 0 ;       % Number of Time series/dimension of shocks
      T = 0 ;       % Number of time periods  
      d = 0 ;          % Length of vector of reduced-form parameters
      
  %% Estimates of VAR parameters (theta)
       AL;          % AL is the estimated lag polynomial
       A0;          % A0 is the estimated constant
       Sigma;       % Sigma is the estimated covariance matrix
       theta;       % Vector of reduced-form parameters
       etahat;      % estimated residuals/forecast errors 
       
 %% Variance 
       Omega;      % asymptotic covariance matrix for reduced form coefficients
       Omegainv;   % Block inverse of Omega (assumes homoskedasticity)

 %% Estimates of IRF
       C;           % C is the estimated reduced form VMA representation, a long matrix (n  x n(MaxHorizons+1) )
       Ccum;        % C is the estimated reduced form VMA representation, a long matrix (n  x n(MaxHorizons+1) )
 
       
 %% todo : check if these properties are necessary
 %% Auxilary properties 
       G;          % Derivatives of C wrt AL.
       Gcum;       %  Derivatives of Ccum wrt AL.
       Vaux;       % auxiliary matrix Vaux such that vec(Sigma)=Vaux vech(Sigma);
       Vaux2;      % auxiliary matrix Vaux2 such that vech(Sigma)=Vaux2 vec(Sigma);

       bootCSigma; % Bootstrap after boostrap samples of of C wrt AL.

         
    end
    
    
    
    
    methods
        function obj = VecAR(label,nLags)
         % This constructor opens a folder ./label/ and reads data
         % and identifying restrictions    
        disp(['Initializing VAR model for ' label]);
            if nargin<2
                obj.nLags = 12 ; 
            else
               obj.nLags = nLags ; 
            end
            disp(['  VAR model has ' num2str(obj.nLags) ' lags.']);
            obj = readTS(obj,label);   % read data   
            % compute LS estimates VAR
            [obj.AL,obj.Sigma,obj.etahat,obj.X,obj.Y,obj.A0] = LS_VAR(obj.data,obj.nLags); 
            [obj.theta,obj.d] = computeTheta(obj.AL,obj.Sigma,obj.nLags,obj.n);
            disp('  LS estimates of VAR model are computed successfully.');
            %% auxiliary matrix Vaux such that: vec(Sigma)=Vaux vech(Sigma)
            [obj.Vaux] = auxiliarymatrix(obj.n);
            %% MA coefficients
            [obj.C,obj.Ccum] = MARep(obj.AL,obj.nLags,obj.MaxHorizons);  
            %% Asymptotic covariance matrix for reduced form coefficients
            [obj.Omega,obj.Omegainv] = CovAhat_Sigmahat(obj.nLags,obj.X,obj.etahat,obj.Vaux,obj.hetscedOmega); 
            %% G matrix: derivative of vec(C) wrt vec(AL) 
            [obj.G,obj.Gcum] = Gmatrices(obj.AL,obj.C,obj.nLags,obj.MaxHorizons,obj.n);  
            
        disp('Initialization of the VAR model is done.');
       
        
        end
    
        function stationarityTest(obj)
% -------------------------------------------------------------------------
% Checks whether reduced-form VAR model is covariance stationary  
% -------------------------------------------------------------------------
            stationarity(obj.AL,obj.nLags,obj.n)  
        end        
        function obj = bab(obj)
% -------------------------------------------------------------------------
% bootstrap after bootstrap   
% See Kilian (1998) ReStat
% 
% Inputs:
% - obj.data  = series: matrix of dimension (T times n) containing the time series
% - obj.nLags = p: model lag order
% - obj.MaxHorizons = hori: model horizon
% - obj.babn1 = n1: number of bootstrap samples before bias correction
% - obj.babn2 = n2: number of bootstrap samples after  bias correction
% Outputs:
% - output: structure containing all reduced-form model objects 
%           after bootstrapping
% 
% This version: March 31, 2015
% Please, cite Gafarov, B. and Montiel-Olea, J.L. (2015) 
% "ON THE MAXIMUM AND MINIMUM RESPONSE TO AN IMPULSE IN SVARS"
% -------------------------------------------------------------------------


%%  read inputs from the VecAR object
series = obj.data;
p      = obj.nLags;
hori   = obj.MaxHorizons ;
n1 = obj.babn1;
n2 = obj.babn2;  
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
    [ALar(:,:,iMC),sigmaAr(:,:,iMC),~] = LS_VAR(seriesMC,p);
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
  [ALar(:,:,iMC),sigmaAr(:,:,iMC),~] = LS_VAR(seriesMC,p);
  [CAr(:,:,iMC)]=MARep(ALar(:,:,iMC),p,hori);
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
series= obj.data;
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
    [~,~,eta,~] = LS_VAR(series_,px); %Apply the function RForm_VAR.m to estimate reduced form parameters 
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
        
        
    end
    
end


function obj = readTS(obj,label) % used in  VecAR constructor
%% this function reads data and labels from a data.csv file in the label folder

%  todo: make filepath an argument 
        %% read data
            obj.data  = csvread([label filesep 'data.csv'],1);
            [obj.T,obj.n] = size( obj.data);                
        %% read names of time series in a cell array
            fid       = fopen([label filesep 'data.csv']);
            obj.names = textscan(fid,[repmat('%[^,],',1,obj.n-1) '%[^,\r\n]'], 1);
            fclose(fid);
            disp('  Data read is succefull.')   ;  % todo: handle exceptions  
            
end

%% functions from folder funcRForm

function [AL,Sigma,eta,X,Y,A0] = LS_VAR(TSL,p) 
% -------------------------------------------------------------------------
% LS_VAR provides  OLS estimators of a VAR(p) model  
% The estimation always includes a constant
% 
% Former name RForm_VAR
%
%
% Inputs:
% - TSL: matrix of dimension (T times n) containing the time series
% - p: number of lags in the VAR model
% Outputs:
% - AL: VAR model coefficients
% - Sigma: covariance matrix of VAR model residuals
% - eta: VAR model residuals
% - X: VAR model regressors
%
% This version: May 4th, 2015
% Last edited by José-Luis Montiel-Olea
% 
% -------------------------------------------------------------------------


%% Definitions
aux = lagmatrix(TSL,1:1:p);
Y = TSL((p+1):end,:); %The rows of this matrix are Y_t'
X = [ones(size(Y,1),1),aux((p+1):end,:)]; clear aux   %The rows of this matrix are [1,X_{t}'] in p.6
 

%% Generate the vec(A(L)) estimators
slopeparameters = (Y'*X)*((X'*X)^(-1)); %contains the nx1 constant vector and AL
AL = slopeparameters(:,2:end); %n x np
A0 = slopeparameters(:,1);     %n x 1
 

%% Covariance matrix
eta = Y'-slopeparameters*(X'); 
Sigma = (eta*eta')/(size(eta,2));
 

end
function [Vaux] = auxiliarymatrix(n)
% -------------------------------------------------------------------------
% computes the auxiliary matrix Vaux such that: vec(Sigma)=Vaux vech(Sigma)
% 
% Inputs:
% - n: number of variables
% Outputs:
% - Vaux: auxiliary matrix
% 
% This version: February 24, 2015
% Please, cite Gafarov, B and Montiel-Olea, J.L. (2015) 
% "ON THE MAXIMUM AND MINIMUM RESPONSE TO AN IMPULSE IN SVARS"
% -------------------------------------------------------------------------

Aux = eye(n);
last = zeros(n^2,1); last(n^2,1)=1;
V(1).V = last;
for j=2:n
    A = eye(j);
    if j<n
        B = [zeros( n*(n-j)+n-j, j);eye(j)];
        for m=2:j
            B = [B; kron(Aux(:,n-j+1),A(m,:))];
        end
        clear A;
        V(j).V = [B,V(j-1).V];
        clear B
    else
        B = eye(j);
        for m=2:j
            B = [B;kron(Aux(:,n-j+1),A(m,:))];
        end
        clear A;
        V(j).V = [B,V(j-1).V];
        clear B
    end    
end
Vaux = V(n).V;
disp('  Vaux is computed: vec(Sigma)=Vaux vech(Sigma).');

end
function [C,Ccum] = MARep(AL,p,hori)
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
% Please, cite Gafarov, B. and Montiel-Olea, J.L. (2015) 
% "ON THE MAXIMUM AND MINIMUM RESPONSE TO AN IMPULSE IN SVARS"
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

disp('  VMA coefficients are computed successfully.')

end
function [OmegaHat,OmegaHatinv] = CovAhat_Sigmahat(p,X,eta,Vaux,hetscedOmega)
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
% This version: March 9, 2015
% Last edited by Bulat Gafarov 
% 

% -------------------------------------------------------------------------


%% suppress default warning for near singularity
warning('off','MATLAB:nearlySingularMatrix')

 

switch hetscedOmega
    case 0

        %% Definitions
        n = size(eta,1);
        XSVARp = X(:,2:end); 
%         XSVARp = X; 
        T1aux = size(eta,2); %This is the number of time periods
        Sigmau = (eta*eta')/(size(eta,2));
        Gamma = XSVARp'*XSVARp/T1aux;
        DK = Vaux; % duplication matrix such that vec(Sigma)=D*vech(Sigma)
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
        
        
    case 1

        %% Definitions
        n = size(eta,1);
        XSVARp = X; 
        matagg = [XSVARp,eta']'; %The columns of this vector are (1;X_t; eta_t)
        T1aux = size(eta,2); %This is the number of time periods
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


        %% Construct the selector matrix Vaux that gives: vech(Sigma)=Vaux*vec(Sigma)
        I = eye(n);
        V = kron(I(1,:),I);
        for i=2:n
            V = [V; kron(I(i,:),I(i:end,:))];
        end


        %% Check bad conditioning
        test = n^2*p + n*(n+1)/2 < T1aux;
        if false(test)
            error('Warning: The covariance matrix is not invertible because n^2p+n(n+1)/2<T!')
        end


        %% This is the estimator for matrix W in Assumption 1
        Mhat = [kron([zeros(n*p,1),eye(n*p)]*(XSVARp'*XSVARp./T1aux)^(-1),eye(n)), zeros((n^2)*p,n^2); ...
                zeros(n*(n+1)/2,((n^2)*p)+n), V];
        OmegaHat = (Mhat)*(WhatAss1)*(Mhat');
        OmegaHatinv = OmegaHat\eye(size(OmegaHat)); 

end


%% Warning for near-singularity
test = rcond(OmegaHat);
if test<eps
    fprintf('... warning: OmegaHat covariance matrix is close to singular!\n')
end
warning('on')
    disp('  Asym. covariance for (AL,Sigma) is computed successfully.')
end
function [G,Gcum] = Gmatrices(AL,C,p,hori,n)
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
% Please, cite Gafarov, B. and Montiel-Olea, J.L. (2015) 
% "ON THE MAXIMUM AND MINIMUM RESPONSE TO AN IMPULSE IN SVARS"
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

disp('  Derivatives of VMA coefficients w.r.t vec(A) are computed successfully.')

end
function simVAR = simulateVAR(AL,T,nMC,etahat)
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
% Please, cite Gafarov, B and Montiel-Olea, J.L. (2015) 
% "ON THE MAXIMUM AND MINIMUM RESPONSE TO AN IMPULSE IN SVARS"
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
% Please, cite Gafarov, B. and Montiel-Olea, J.L. (2015) 
% "ON THE MAXIMUM AND MINIMUM RESPONSE TO AN IMPULSE IN SVARS"
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
function [theta,d] = computeTheta(AL,Sigma,p,n)

%% this funciton computes the vector of the reduced form parameters, theta vector

Ident = eye(n);
Vaux2 = kron(Ident(1,:),Ident);
for i=2:n
    Vaux2 = [Vaux2; kron(Ident(i,:),Ident(i:end,:))];
end
vecAL = reshape(AL,[(n^2)*p,1]);
vechSigma = Vaux2*reshape(Sigma,[n^2,1]);
theta = [vecAL; vechSigma]; 
d = length(theta);
end


%********************************************************
%********************************************************
%********************************************************
%********************************************************