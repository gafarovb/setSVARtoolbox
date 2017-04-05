classdef SVAR < estimatedVecAR 
    
    %%        SVAR class describes a set-identified SVAR model and offers tools
    %       to construct point estimates and conduct inteference on the impulse
    %       response functions (IRF).
    %
    %       This version: 3.9.2017
    %       Last edited by Bulat Gafarov
    %%

    
    properties (Access = public)
        label = 'Unknown' ; % Model label, e.g. MSG 
        spec ;      % intermediate computations based restrictions.
        ID   =[] ; %% An object with restrictions

    end
    
    properties (Access = private)
        optimiaztionProblems;
        % todo: make reduced VAR a  property!
    end
        
    methods
        function obj = SVAR
            obj@estimatedVecAR; % create a reduced form VAR model
            obj.label = obj.config.label;
            obj.ID = IDassumptions( obj.config.assumptionsFilename); % read ID assumptions from a file
            obj.optimiaztionProblems = optimizationProblems( obj);  
         end
        
        %%
        function maxIRFcollection = onesidedUpperIRFHatAnalytic(obj)
            maxBoundsMatrix = obj.optimiaztionProblems.getMaxBounds;    
            maxIRFcollection = IRFcollection(maxBoundsMatrix, obj.getNames, 'max point estimates') ;
        end
        function minIRFcollection = onesidedLowerIRFHatAnalytic(obj)
            minBoundsMatrix = obj.optimiaztionProblems.getMinBounds;    
            minIRFcollection = IRFcollection(minBoundsMatrix, obj.getNames, 'min point estimates') ;
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
        function cumOrNot = cum(obj)
            cumOrNot = obj.config.cum;
        end
        function restMat = getRestMat(obj)
           restMat = obj.ID.getRestMat;
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
vecFromVech  = SVARobj.vecFromVech;
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
            %% TODO
            H = Hin(:,m,t) ; %%% THIS IS WRONG.
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
                (B/2)*(kron(H'*sigmaInv, H'*sigmaInv))*vecFromVech];
            
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
restMat = SVARobj.ID.getRestMat;
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
