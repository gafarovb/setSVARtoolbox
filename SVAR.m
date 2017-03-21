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
        spec;      % intermediate computations based restrictions.

    end
    
    properties (Access = private)
        ID   =[] ; %% A class with restrictions
        SR  ; % matrix with sign restrictions
        ZR  ; % matrix with zero restrictions
        GSR ; % Derivatives of sign restricted IRF
        GZR ; % Derivatives of zero restricted IRF
    end
    
    %% *****************************************************************
    %  *****************************************************************
    %% *****************************************************************
    
    methods
        function obj = SVAR(label,nLags)
            %  INPUTS:
            %  label is the name of the model
            %  nLags is the number of lags in the model
            obj@VecAR(label,nLags); % create a reduced form VAR model
            obj.label = label;
            obj.ID = restrictions(label); % read restrictions from a file
            obj = obj.separateSandZ;             % creates a specificaiton structure to characterize the restrictions
        end
        
        %%
        function [IRFs,solution] = onesidedIRFHatAnalytic(obj,minmax)
            % todo: refactor this method
            
            
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
                    maskCol = obj.ID.maskColMin;
                    fiq = spec.FiqMin;
                    sortDirection = 'ascend';
                case 'max'
                    minInd = 1;
                    maskCol = obj.ID.maskColMax;
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
            IRFs = IRFcollection(value,obj.names,[minmax, ' point estimates']) ;
            solution = struct( ...
                'Value',value, ...   % lower or upper bounds for every TS/horizon
                'arg',arg, ...       % argmax vectors for every TS/horizon
                'Value2',Bound2,...  % second largest values that don't violate the sign constraints and the "nestedness" rule
                'HAr',HArnew,...
                'valAr',valArnew);
            
        end
        
        
        function obj = separateSandZ(obj)
            %% -------------------------------------------------------------------------
            % This function constructs separate matrices with sign and zero restrictions
            %
            % Inputs:
            % - obj: SVAR object
            %
            % Outputs:
            % - specification: structure containing all information from the
            %                  restriction matrix
            %
            % This version: March 21th, 2017
            % Last edited :  Bulat Gafarov
            % -------------------------------------------------------------------------
            
            
            
            
            %% interface
            restMat = obj.ID.restMat;
            nSignRestrictions = countSignRestrictions(obj.ID);
            nZeroRestrictions = countZeroRestrictions(obj.ID);
            obj.ID = constructSelectors( obj.ID,obj.n);
            n = obj.n;
            [nG,mG,MaxHorizons] = size(obj.G);
            VMA    = reshape([eye(n),obj.C]   ,[n,n,MaxHorizons]);
            VMAcum = reshape([eye(n),obj.Ccum],[n,n,MaxHorizons]);
            
            %% allocate memory
            obj.SR = zeros(nSignRestrictions,n); % matrix with sign restrictions
            obj.ZR = zeros(nZeroRestrictions,n); % matrix with zero restrictions
            obj.GSR = zeros(nG,mG,nSignRestrictions); % Derivatives of sign restricted IRF
            obj.GZR = zeros(nG,mG,nZeroRestrictions); % Derivatives of zero restricted IRF
            
            iSR = 1;
            iZR = 1;
            
            for i=1:(nSignRestrictions+nZeroRestrictions) % loop through all restricions
                
                currentRestriction = (1-restMat(i,4)) *     VMA(restMat(i,1),:,restMat(i,2)+1) +...
                    restMat(i,4)  *  VMAcum(restMat(i,1),:,restMat(i,2)+1);
                currentRestrictionDerivative = (1-restMat(i,4)) *    obj.G(:,:,restMat(i,2)+1) +...
                    restMat(i,4)  * obj.Gcum(:,:,restMat(i,2)+1);
                
                if  restMat(i,3)==0  % check if it is a zero restrictions
                    obj.ZR(iZR,:)    = currentRestriction;
                    obj.GZR(:,:,iZR) = currentRestrictionDerivative;
                    
                    iZR=iZR+1;
                else    % sign constraints
                    
                    obj.SR(iSR,:) = currentRestriction;
                    obj.GSR(:,:,iSR) = currentRestrictionDerivative;
                    %% Convert restrictions to cannonical form
                    obj.SR(iSR,:)     = restMat(i,3)* obj.SR(iSR,:) ;
                    obj.GSR(:,:,iSR)  = restMat(i,3)* obj.GSR(:,:,iSR) ;
                    
                    iSR=iSR+1;
                end
            end
            %% legacy interface
            
            obj.spec = struct(...
                'Niq',nSignRestrictions,...    % number of sign (inequality) restrictions
                'Neq',nZeroRestrictions,...    % number of zero (equality) restrictions
                'eiq',obj.ID.selectorSR,...    % 'e' vectors of inequality restrictions
                'eeq',obj.ID.selectorZR,...    % 'e' vectors of equality restrictions
                'Giq',obj.GSR,...    % G matrix associated with inequality constraints
                'Geq',obj.GZR,...    % G matrix associated with   equality constraints
                'Ziq',(obj.SR)',...   % sign restrictions
                'Zeq',(obj.ZR)',...   % zero restrictions
                'maskColMax',obj.ID.maskColMax,...    % IRF to skip
                'maskColMin',obj.ID.maskColMin,...    % IRF to skip
                'FiqMax', obj.ID.negativeSR,...  % matrix with indexed negative sign restrictions
                'FiqMin', obj.ID.positiveSR,...   % matrix with indexed positive sign restrictions
                'restMat',obj.ID.restMat);   % input restriction matrix
        end
        function IRFcol     = onesidedIRFCSAnalytic(obj,minmax,level)
            
            %% baseline analytical algorithm in GM 2014
            
            [~, pointEstimates]  = onesidedIRFHatAnalytic(obj,minmax);
            switch(minmax)
                case 'min'
                    IRFmatrix = pointEstimates.Value - stdIRF(obj,pointEstimates) * norminv(level,0,1);
                case 'max'
                    IRFmatrix = pointEstimates.Value + stdIRF(obj,pointEstimates) * norminv(level,0,1);
                otherwise
                    error('Choose min or max');
            end;
            IRFmatrix = enforceIRFRestrictions(obj,IRFmatrix) ;
            IRFcol    = IRFcollection(IRFmatrix,obj.names,[minmax,' analytic CS with p=',num2str(level)]);
        end   
        function objSimulated = resampleTheta(obj,seedMC)
            objSimulated = resampleTheta@VecAR(obj,seedMC);
            objSimulated = objSimulated.separateSandZ; 
        end
        function n = nTS(obj) % fixme
            % Funciton nTS returns the number of time series
            n = obj.n;
        end
        function cumOrNot = cum(obj)
            cumOrNot = obj.config.cum;
        end
        
    end
    
end

%--------------- Folder ./funcSForm -----------------------


function [CV] = stdIRF(SVARobj,solution)


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


%% compute std errors using Delta method
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

function IRF = enforceIRFRestrictions(SVARobj,IRF)
%% This function enforces sign and zero restrictions specified in SVARobj
%  on IRF
%
% last modified : March 15 2017
% by Bulat Gafarov

%% read variables from the input objects
nShocks = SVARobj.n;
restMat = SVARobj.ID.restMat;
nRestrictions = size(restMat,1);
MaxHorizons = SVARobj.config.MaxHorizons;

%% allocate memory
restrictedBelow = false(nShocks,MaxHorizons+1);
restrictedAbove = false(nShocks,MaxHorizons+1);

for i = 1:nRestrictions
    iVariable        = restMat(i,1);
    iHorizon         = restMat(i,2)+1;
    iRestrictionType = restMat(i,3);
    switch iRestrictionType
        case 0
            restrictedBelow(iVariable,iHorizon) = true;
            restrictedAbove(iVariable,iHorizon) = true;
        case 1
            restrictedBelow(iVariable,iHorizon) = true;
        case -1
            restrictedAbove(iVariable,iHorizon) = true;
        otherwise
            ME = MException('enforceIRFRestrictions:incorrectRestriction', ...
                'Unknown restiction type: %d',iRestrictionType);
            throw(ME)
    end
end

IRF(restrictedAbove) = min(0,IRF(restrictedAbove));
IRF(restrictedBelow) = max(0,IRF(restrictedBelow));

end
