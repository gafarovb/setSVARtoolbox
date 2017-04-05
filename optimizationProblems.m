classdef optimizationProblems < handle
    %OPTIMIZATIONPROBLEMS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        subproblems;
        
        objectiveFunctions;
        
        Sigma;
        
        linearConstraints;
        nInequalities;
        nEqualities;
        
        config;
    end
 
    methods
        function obj = optimizationProblems( objSVAR)
            obj.nInequalities = objSVAR.ID.countSignRestrictions;  % number of sign restrictions
            obj.nEqualities   = objSVAR.ID.countZeroRestrictions;  % number of zero restrictions
            obj.config = objSVAR.getConfig;
            obj.Sigma   = objSVAR.getSigma;
            obj.objectiveFunctions = objSVAR.getReducedFormIRF;
            obj.linearConstraints =  objSVAR.ID.generateLinearConstraints(objSVAR); 
            
            obj.subproblems = initializeSubproblems(obj);
            
        end
        function SigmaSqrt = getSQRTSigma(obj)
            SigmaSqrt = (obj.Sigma)^(1/2);    %Sigma^(1/2) is the symmetric sqrt of Sigma.
        end
        function subProblems = initializeSubproblems(obj)
            degreesOfFreedom = obj.getN-obj.nEqualities-1;
            MaxNumberOfActiveSR = min( obj.nInequalities, degreesOfFreedom);
            % SR stands for Sign Restrictions 
            
            iSubProblem = 1; 
            activeSet = false(obj.nInequalities,1);
            subProblems = subproblemsGivenActiveSet(activeSet,obj);
            
            for nActiveSR = 1:MaxNumberOfActiveSR
                
                indicesOfActiveSRForSubsetsOf_nActive_elements = ...
                    obj.getAllsubsetsSize_k_from_N(  obj.nInequalities,  nActiveSR);
                nSubProblemsOfSize_nActiveSR = size( indicesOfActiveSRForSubsetsOf_nActive_elements ,1);
                
                for ix = 1 : nSubProblemsOfSize_nActiveSR
                    iSubProblem = iSubProblem + 1;
                    activeSet = false( obj.nInequalities,  1);
                    activeSet(indicesOfActiveSRForSubsetsOf_nActive_elements(ix,:)) = true; % create an index mask for a given active set of constraint with index jx

                    subProblems(iSubProblem) = subproblemsGivenActiveSet( activeSet, obj);
                end
            end
            
        end
        function linearConstraints = getLinearConstraints(obj)
            linearConstraints = obj.linearConstraints;
        end
        function objectiveFunctions = getObjectiveFunctions(obj)
            objectiveFunctions = obj.objectiveFunctions;
        end
        function configHandle = getConfig (obj)
            configHandle = obj.config;
        end
        function nSubProblems = countSubproblems(obj)
            nSubProblems = size( obj.subproblems,2);
        end
        function maxHorizons = getHorizons(obj)
            maxHorizons = size(obj.objectiveFunctions,3);
        end
        function nShocks = getN(obj)
            nShocks = size(obj.objectiveFunctions,1);         
        end
        function maxBounds = getMaxBounds(obj)
            maxBoundsForSubproblems  =  zeros(obj.getN,obj.getHorizons,obj.countSubproblems);
            nSubProblems = obj.countSubproblems;
            
            for i = 1:nSubProblems
                maxBoundsForSubproblems(:,:,i) = obj.subproblems(i).maxBounds;
            end
            maxBounds = max(maxBoundsForSubproblems,[],3);
        end
        
        function minBounds = getMinBounds(obj)
            minBoundsForSubproblems  = zeros(obj.getN,obj.getHorizons,obj.countSubproblems);
            nSubProblems = obj.countSubproblems;
            
            for i = 1:nSubProblems
                minBoundsForSubproblems(:,:,i) = obj.subproblems(i).maxBounds;
            end
            minBounds = min( minBoundsForSubproblems,[],3);
        end
    end
    
    methods (Static)
        function indicesOfActiveSRForSubsets = getAllsubsetsSize_k_from_N(N,k)
            if k == 0  % no active restrictions
                indicesOfActiveSRForSubsets = []; 
            else
                indicesOfActiveSRForSubsets = combnk(1: N, k); 
            end
            
        end
    end
end

