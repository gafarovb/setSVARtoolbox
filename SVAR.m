classdef SVAR < VecAR
    
%%        SVAR class describes a set-identified SVAR model and offers tools
%       to construct point estimates and conduct inteference on the impulse 
%       response functions (IRF).
%
%       This version: 3.9.2017
%       Last edited by Bulat Gafarov 
%%    
    properties (Access = public) 
      label = 'Unknown' ; %% Model label, e.g. MSG 
      ID   =[] ; %% A class with restrictions
    end
    
    properties (Access = private, Constant)      % todo: make it a cofiguration file
        %% Specify shock
        % solve for IRF bounds of specified shock 
        % (by default, compute bounds for all shocks under gridsearch algorithm)
    shock = 1; 
%% further model options
    coverage = 0; % 0=off, 1=on: compute MC coverage frequency
    cluster = 0;  % 0=off, 1=on: use cluster to compute MC coverage frequency
    caliProj = 0; % 0=off, 1=on: 'fake' calibration of Projection CS, 2=on: calibration using IRF rotations 
    sobol = 0;    % 0=off, 1=on: use Sobol sequences instead of random numbers
    bootstrap = 0; % 0=draw reduced-form parameters, 1=bootstrap data and re-estimate reduced-form parameters (instead of drawing them directly)
    dif_calibration = 0;
    nd = 1e5;     % number of (accepted) grid search draws
    nBoot = 1000; % Number of bootstrap/projeciton draws
    nMC = 2000;   % Number of MC draws for coverage frequency
    nblocks = 50; % Number of workers on cluster
    cum = '';  % 'cum' = cumulative IRF responses 
    level = 0.68; % confidence level of bounds on bounds 
    mineig = 0.001; % lower bound on smallest eigenvalue of sigma
    maxmod = 0.990; % upper bound on largest modulus of eigenvalue of companion matrix of A
 
 
   end
    
    properties (Access = private)
        spec; % restrictions specifications.
    end
    
%% *****************************************************************  
%  *****************************************************************  
%  *****************************************************************  
%  *****************************************************************      
%% *****************************************************************  

    methods
      function obj = SVAR(label,nLags)

         %  INPUTS: 
         %  label is the name of the model
         %  nLags is the number of lags in the model
      
         obj@VecAR(label,nLags); % create a reduced form VAR model
         obj.label = label ;
         obj.ID = restrictions(label); % read restrictions from a file 
         obj = cnstr(obj);             % creates a specificaiton structure to characterize the restrictions
      end
      
 %% 
 
      function solution = onesidedIRFHatAnalytic(obj,minmax)
     
% -------------------------------------------------------------------------
% Computes point estimates for lower/upper bounds on the structual IRFs 
% with m<=n-1 zero-restrictions and sign-restrictions. The sign restriction 
% need not be contemporaneous.
% Hint: One can generate lower bounds by supplying -Z instead of Z and by
% multiplying the output by -1.
% One can compare BoundU and BoundU2 to decide if a given value function
% can potentially have at least two argmaxima. In that case the
% corresponding components will be less than some treshold value epsilon2
% This comparison only make sense if there are 2 or more sign restrictions
% in the model.
%
% Former name : BoundsGM(minmax,RVAR,spec,opt)
% 
% Inputs:
% - obj: SVAR object
% - minmax: 'max' ('min') for upper (lower) bound on structural IRFs
% Outputs:
% - solution: structure with bounds on structural IRFs
%  
% This version: March 10, 2017
% Please, cite Gafarov, B. and Montiel-Olea, J.L. (2015) 
% "ON THE MAXIMUM AND MINIMUM RESPONSE TO AN IMPULSE IN SVARS"
% -------------------------------------------------------------------------

spec = obj.spec;
%% Read input structure
% Min or max bounds
switch(minmax)    
    case 'min'
        minInd = -1;
        maskCol = spec.maskColMin;
        fiq = spec.FiqMin; 
        sortDirection = 'ascend'; 
    case 'max'
        minInd = 1;
        maskCol = spec.maskColMax;
        fiq = spec.FiqMax;
        sortDirection = 'descend';
    otherwise
        error('Choose min or max');
end;

%% Simple or cumulative IRF
switch(obj.cum) 
    case 'cum'
        C = obj.Ccum;
    otherwise
        C = obj.C;
end;


Sigma = obj.Sigma;
Zeq   = spec.Zeq';
Z     = spec.Ziq';


%% Definitons
n = size(C,1);                %Gives the dimension of the VAR
hori = size(C,2)/n;           %Gives the horizon of IRFs
largeNum = 100000*minInd;   % fixme 
C = [eye(n),C]; 
Chataux = reshape(C,[n,n,(hori+1)]);
sigmaSqrt = (Sigma)^(1/2);    %Sigma^(1/2) is the symmetric sqrt of Sigma.


%% Allocate memory 
Bound2  = zeros(n,hori+1);
[Ms,~] = size(Z); % number of sign restrictions
[Mz,~] = size(Zeq); % number of zero restrictions



% -------------------------------------------------------------------------
%% Define nested functions 

function [Boundloc,BoundlocNeg,Hloc,valloc,conviol,conviolNeg]=restricted(active,fiq) 
% this subfunction computes the upper bounds BoundU and the corresponding
% argmaxes HU given that Zactive constraints are active and sign constraints Z are imposed    


    %% initialize local variables
    Zactive  = [Zeq; Z(active,:)];
    Znact    = Z(~active,:);
    Boundloc    = zeros(n,hori+1);
    BoundlocNeg = zeros(n,hori+1);
    conviol    = zeros(n,hori+1);
    conviolNeg = zeros(n,hori+1);
    valloc    = zeros(n,hori+1);
    Hloc = zeros(n,n,hori+1);  

    
    %% 'effective' covariance matrix under active constraints
    %  Compare with Lemma 1
    if ~isempty(Zactive) % if there are active zero constraints 
        M = eye(n) - ( (sigmaSqrt * Zactive')*((Zactive * Sigma * Zactive')^(-1))* (Zactive * sigmaSqrt)) ; % compare with  Proposition 1
        sigmaM = sigmaSqrt * M * sigmaSqrt;   % "Effective" covariance matrix given the set Zactive
    else % if there are no active zero constraints, follow the formula for the  unresrticted case
       sigmaM = Sigma ;
    end
    
    
    for ix=1:hori+1  % loop over IRF horizon
        
        %% bounds under active constraints
        Htmp = sigmaM*Chataux(:,:,ix)';
        % Compute max value under active constraints
        % - abs is used to avoid complex numbers with Imaginary part equal to Numerical zero
        boundtmp = abs((diag(Chataux(:,:,ix) * Htmp)).^.5);  
        % Compute argmax under active constraints
        %  - bsxfun divides by zeros without warning
        Hloc(:,:,ix) = bsxfun(@rdivide,Htmp,boundtmp'); 
%         Hloc(:,:,ix) = Htmp./(ones(n,1)*boundtmp'); 

        
%         Hloc(:,:,ix)=spec.Bounds_ex3.mini.arg(:,:,2);
        
        
        %% Check whether signs restrictions are satisfied. Compare with Proposition 1
        mask  = true(n,1)  ;
        slackness = zeros(size(Znact,1),n); 
        for iFiq=1:size(fiq,1)
           if   (fiq(iFiq,2)+1)==ix
               if active( fiq(iFiq,4) )
                   mask(fiq(iFiq,1))=false;
               else
                   columnWithOne = zeros(size(Z,1),1);
                   columnWithOne(fiq(iFiq,4))=abs(largeNum);
                   columnWithOne = columnWithOne(~active,:); 
                   slackness(:,fiq(iFiq,1))= columnWithOne ; %this matrix relaxes 'self'- sign restrictions 
               end
           end 
        end
        if isempty(Znact)  % if there are no sign restrictions, the idicator functions are equal to 1
            Boundloc(:,ix)    =   boundtmp ; 
            BoundlocNeg(:,ix) = - boundtmp ; 
        else  % iff  some sign restriction is violated by a margin of a numerical error of 1.e-6
              % the corresponding indicator function is equal to 0  
            Boundloc(:,ix)    =   boundtmp - largeNum * (ones(n,1) - ( all(Znact*Hloc(:,:,ix)+slackness>-1.0e-6,1)' & mask ) );
            BoundlocNeg(:,ix) = - boundtmp - largeNum * (ones(n,1) - ( all(Znact*Hloc(:,:,ix)-slackness< 1.0e-6,1)' & mask ) );
        end
        conviol(:,ix)    = min( Znact*Hloc(:,:,ix))'<-1e-6;
        conviolNeg(:,ix) = min(-Znact*Hloc(:,:,ix))'<-1e-6;
        valloc(:,ix)    =   boundtmp;
    end
end
% -------------------------------------------------------------------------


%% This segment implements the formula from Proposition 1 and 2


%% compute number of combinations of active sign restrictions
nComb = 1;
for iBin = 1:min(n-1,Ms)
    nComb = nComb + nchoosek(Ms,iBin);
end


%% initialize
BoundsAr = - largeNum*ones(n,hori+1,nComb*2);    
HAr      = zeros(n,n,hori+1,(nComb)*2);  
Boundlx  = zeros(n,hori+1);
valAr = zeros(n,hori+1,nComb*2);    
ConViol = zeros(n,hori+1,nComb*2);    


%% loop through all possible cardinalities of subsets of active sign restrictions 
lg = -1; % total index for all cardinalities of active sets 
for lx = 0:min(Ms,n-Mz-1)  
    
    
    %% count number of subsets of given length lx
    if lx==0  % no active restrictions
        combinations = []; % set of combinations is empty!
        lN = 1; % we need to consider only the unrestricted case
    else
        combinations = combnk(1:Ms,lx); % matrix with indexes of all subsets of length lx
        [lN,~] = size(combinations);   % number of combinations of length lx to consider 
    end
    
    
    %% loop through all possible subsets of sign restrictions of lengths lx
    ll = lg+2; % the first index of active set of length lx
    for jx = 1:lN 
        activeSet = false(Ms,1);  
        if lx>0  % active set could be empty
            activeSet(combinations(jx,:)) = true; % create an index mask for a given active set of constraint with index jx
        end
        lg = lg+2;   
        % activate only restrictions from activeSet 
%         [BoundsAr(:,:,lg),BoundsAr(:,:,lg+1),HAr(:,:,:,lg)] = restricted(activeSet,fiq);  
        [BoundsAr(:,:,lg),BoundsAr(:,:,lg+1),HAr(:,:,:,lg),valAr(:,:,lg),ConViol(:,:,lg),ConViol(:,:,lg+1)] = restricted(activeSet,fiq);  
        HAr(:,:,:,lg+1) = -HAr(:,:,:,lg);
        valAr(:,:,lg+1) = -valAr(:,:,lg);
    end
    lu = lg+1; % the last index of active set of length lx

    
    %% check for multiplicity of solutions, see Proposiion 2
    [BoundlxSort,~] = sort(BoundsAr(:,:,ll:lu),3,sortDirection); 
    BoundlxMAX = (lx+1)*(BoundlxSort(:,:,1)>0); % check if we found max that doesn't violate the sign constraing for the cardinality lx  
    if lu>ll  %if there are different active set of cardinality lx, ...
        BoundlxMAX2 = BoundlxSort(:,:,2); %  store second largest value, possibly zero
    else
        BoundlxMAX2 = zeros(n,hori+1); %  otherwise store zero
    end
    % check for "nestedness" of the second largest value
    Boundlx(Boundlx==0)     = BoundlxMAX(Boundlx==0); % this array stores cardinalities for the maximums that we have found already
    Bound2(Boundlx==(lx+1)) = BoundlxMAX2(Boundlx==(lx+1)); % this operation stores  2nd largest non-zero for cardinalities lx 
                                                            % only if we already have found a corresponding maximum of the same cardinality   
end


%% generate output
% Find maximum over all combinations of active sets
[BoundsArSort,indexSet] = sort(BoundsAr,3,sortDirection); 
value = BoundsArSort(:,:,1);  % maximum bounds
arg =  zeros(n,n,hori+1);     % argmax vectors for every TS/horizon
for k = 1:(hori+1) 
    for ts=1:n
        arg(:,ts,k) = HAr(:,ts,k,indexSet(ts,k,1));
    end
end  
value(maskCol) = ((minInd*value(maskCol))<0).*value(maskCol);

valArlong = reshape(valAr,[n,hori+1,2,nComb]);
[~,ind] = max(minInd*valArlong,[],3);
valArnew = zeros(n,hori+1,nComb);
HArlong = reshape(HAr,[n,n,hori+1,2,nComb]);
HArnew = zeros(n,n,hori+1,nComb);
for i=1:n
    for h=1:hori+1
        for j=1:nComb
            valArnew(i,h,j) = valArlong(i,h,ind(i,h,1,j),j);
            HArnew(:,i,h,j) = HArlong(:,i,h,ind(i,h,1,j),j);
        end
    end
end

solution = struct( ...
           'Value',value, ...   % lower or upper bounds for every TS/horizon
           'arg',arg, ...       % argmax vectors for every TS/horizon
           'Value2',Bound2,...  % second largest values that don't violate the sign constraints and the "nestedness" rule
           'HAr',HArnew,...    
           'valAr',valArnew);    

      end

      function solution = onesidedIRFCSAnalytic(obj,minmax,level)
 
%% baseline analytical algorithm in GM 2014
     solutionHat  = onesidedIRFHatAnalytic(obj,minmax);
     dltCVU = std(obj,solutionHat); 
% todo switch here min/max
     solution  = Bounds.Value + dltCVU*norminv(1-(1-level)/2,0,1);
 
  
% todo : make the following code a separate function
%% ensure identifying restrictions satisfied


%     zeroLB = zeros(RVARhat.n,opt.horizon+1);
%     zeroUB = zeros(RVARhat.n,opt.horizon+1);
%     for i=1:size(Bounds.spec.restMat,1)
%         iv = Bounds.spec.restMat(i,1);
%         ih = Bounds.spec.restMat(i,2)+1;
%         ir = Bounds.spec.restMat(i,3);
%         if ir==0
%             zeroLB(iv,ih) = 1;
%             zeroUB(iv,ih) = 1;
%         elseif ir==1
%             zeroLB(iv,ih) = 1;
%         elseif ir==-1
%             zeroUB(iv,ih) = 1;
%         end        
%     end
%     for t=1:opt.horizon+1
%         dltU(zeroUB(:,t)==1,t) = min(0,dltU(zeroUB(:,t)==1,t));
%         dltU(zeroLB(:,t)==1,t) = max(0,dltU(zeroLB(:,t)==1,t));
%         dltL(zeroUB(:,t)==1,t) = min(0,dltL(zeroUB(:,t)==1,t));
%         dltL(zeroLB(:,t)==1,t) = max(0,dltL(zeroLB(:,t)==1,t));
%     end
%     
    
    
      end
      
      function n = nTS(obj) % fixme
         % Funciton nTS returns the number of time series
          n = obj.n;
      end 
      
    end
    
end

%--------------- Folder ./funcSForm -----------------------

function obj = cnstr(obj) 
% -------------------------------------------------------------------------
% This function constructs a structure with the specification of the
% restrictions and auxilary matrices
% 
% Inputs:
% - obj: SVAR object
%
% Outputs:
% - specification: structure containing all information from the
%                  restriction matrix
% 
% This version: May 11th, 2015
% Please, cite Gafarov, B. and Montiel-Olea, J.L. (2015) 
% "ON THE MAXIMUM AND MINIMUM RESPONSE TO AN IMPULSE IN SVARS"
% -------------------------------------------------------------------------

% todo : refactor this code

%% Definitions
restMat = obj.ID.restMat;
Niq = sum(abs(restMat(:,3))); % number of sign restrictions
Neq = sum(restMat(:,3)==0);   % number of zero restrictions
% Nfiq =  sum(abs(restMat(:,3))&restMat(:,2)>0); % number of non-contemporaneous sign restrictions
n = obj.n; 
[nG,mG,hori] = size(obj.G);
C    = reshape([eye(n),obj.C]   ,[n,n,hori]);
Ccum = reshape([eye(n),obj.Ccum],[n,n,hori]);
Fiq  = zeros(Niq,4);% matrix that orders sign restricions
Ziq = zeros(Niq,n); % matrix with sign restrictions
Zeq = zeros(Neq,n); % matrix with zero restrictions
Giq = zeros(nG,mG,Niq); % Derivatives of sign restricted IRF
Geq = zeros(nG,mG,Neq); % Derivatives of zero restricted IRF
eiq = zeros(n,Niq);
eeq = zeros(n,Neq);
Emat = eye(n);

iiq = 1;
ieq = 1;
lf  = 0;
for i=1:(Niq+Neq) % loop through all restricions
    if  restMat(i,3)==0  % check if it is a zero restrictions
        Zeq(ieq,:) = (1-restMat(i,4)) *    C(restMat(i,1),:,restMat(i,2)+1) +...
                        restMat(i,4)  * Ccum(restMat(i,1),:,restMat(i,2)+1);
        Geq(:,:,ieq) = (1-restMat(i,4)) *    obj.G(:,:,restMat(i,2)+1) +...
                           restMat(i,4) * obj.Gcum(:,:,restMat(i,2)+1);      
        eeq(:,ieq) = Emat(:,restMat(i,1));
                     
        ieq=ieq+1;     
    else    % sign constraints
        lf=lf+1;
        Fiq(lf,1:3)=restMat(i,1:3);
        Fiq(lf,4)= iiq;
        Ziq(iiq,:) = (1-restMat(i,4)) *    C(restMat(i,1),:,restMat(i,2)+1) +...
                         restMat(i,4) * Ccum(restMat(i,1),:,restMat(i,2)+1);
        Giq(:,:,iiq) = (1-restMat(i,4)) *    obj.G(:,:,restMat(i,2)+1) +...
                           restMat(i,4) * obj.Gcum(:,:,restMat(i,2)+1);      
        Ziq(iiq,:)    = restMat(i,3)* Ziq(iiq,:) ;
        Giq(:,:,iiq)  = restMat(i,3)* Giq(:,:,iiq) ;
        eiq(:,iiq) = Emat(:,restMat(i,1));
        iiq=iiq+1;     
    end
    
end


FiqMax = Fiq(Fiq(:,3)<0,:); % matrix with indexed negative sign restrictions
FiqMin = Fiq(Fiq(:,3)>0,:); % matrix with indexed positive sign restrictions
if ~isempty(FiqMax)
    maskColMax  = logical( full(sparse(FiqMax(:,1),        1+FiqMax(:,2),1,n,hori+1))); 
else
    maskColMax  = false(n,hori+1)  ;
 
end
if ~isempty(FiqMin)
    maskColMin  = logical( full(sparse(FiqMin(:,1),        1+FiqMin(:,2),1,n,hori+1))); 
    
else
    maskColMin  = false(n,hori+1)  ;  
end


obj.spec = struct(...
       'Niq',Niq,...    % number of sign (inequality) restrictions
       'Neq',Neq,...    % number of zero (equality) restrictions
       'eiq',eiq,...    % 'e' vectors of inequality restrictions
       'eeq',eeq,...    % 'e' vectors of equality restrictions
       'Giq',Giq,...    % G matrix associated with inequality constraints
       'Geq',Geq,...    % G matrix associated with   equality constraints
       'Ziq',Ziq',...   % sign restrictions
       'Zeq',Zeq',...   % zero restrictions
       'maskColMax',maskColMax,...    % IRF to skip
       'maskColMin',maskColMin,...    % IRF to skip
       'FiqMax',FiqMax,...  % matrix with indexed negative sign restrictions
       'FiqMin',FiqMin,...   % matrix with indexed positive sign restrictions
       'restMat',restMat);   % input restriction matrix
end

function [CV] = std(SVARobj,solution)
% -------------------------------------------------------------------------
% Computes the level-% upper confidence band for the max impulse using 
% the Delta-Method. 
% 
% Inputs:
% Outputs
% -------------------------------------------------------------------------


%% Read input structure
cum   = SVARobj.cum;
Bound = solution.Value; % lower or upper bounds for every TS/horizon
Hin   = solution.arg;   % argmax vectors for every TS/horizon
 
valAr = solution.valAr; % 
% Simple or cumulative IRF
switch(cum)    
    case 'cum'
        C = SVARobj.Ccum; % fixme , do not use C and G  
        G = SVARobj.Gcum;
    otherwise
        C = SVARobj.C;
        G = SVARobj.G;
end;
Sigma = SVARobj.Sigma;
Omega = SVARobj.Omega;
Vaux  = SVARobj.Vaux;
T     = SVARobj.T;
d     = SVARobj.d;
spec  = SVARobj.spec;
 


%% Definitons
n    = size(C,1);
hori = size(C,2)/n;
C    = [eye(n),C];
C    = reshape(C,[n,n,(hori+1)]);
e    = eye(n);
sigmaInv = Sigma^(-1);
nComb = size(valAr,3);


%% allocate memory
CV  =  zeros(n,hori+1);


%% compute confidence bands using Delta method
for m=1:n
    for t=1:hori+1 
        
        CVlong = zeros(nComb,1);
        
        for j=1:nComb

            H = Hin(:,m,t) ;
            B = Bound(m,t);
             
            %% compute adjustment for non-contemporaneous restrictions 
            maskActIneq = (spec.Ziq'*H<1e-5);
            Zact = [spec.Ziq(:,maskActIneq)';spec.Zeq'];
            Eact = [spec.eiq(:,maskActIneq)  spec.eeq ];
            Gact = cat(3, spec.Giq(:,:,maskActIneq'), spec.Geq);
            w = (Zact * Sigma *  Zact')\Zact * Sigma *C(:,:,t)'*e(:,m);
            adjustment = 0;
            for iz=1:size(w) 
                adjustment = adjustment - w(iz)*kron(H',Eact(:,iz)')*Gact(:,:,iz);
            end

            %% gradient of parameter vector theta
             Grad = [kron(H',e(:,m)')*G(:,:,t)+adjustment,...
                    (B/2)*(kron(H'*sigmaInv, H'*sigmaInv))*Vaux];

            %% gradient of argmax
             CVlong(j) = abs((Grad*Omega*Grad')^.5)/(T^.5);
            
        end
        
        CV(m,t) = max(CVlong);
        
    end
end

CV(isnan(CV)) = 0 ; % report 0 if an error has happend


end