classdef optimizationProblems < handle
    %OPTIMIZATIONPROBLEMS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        nInequalities;
        nEqualities;
        nShocks;
        subproblems;
        objectiveFunctions;
        linearConstraints;
    end
    
    methods
        function obj = optimizationProblems( objSVAR)
            obj.nShocks = objSVAR.getN;
            obj.nInequalities = objSVAR.ID.countSignRestrictions;  % number of sign restrictions
            obj.nEqualities   = objSVAR.ID.countZeroRestrictions;  % number of zero restrictions
            
            obj.objectiveFunctions = objSVAR.getReducedFormIRF;
            obj.linearConstraints =  objSVAR.ID.generateLinearConstraints(objSVAR); 
            
            obj.subproblems = initializeSubproblems(obj);
            
        end
        
        function subProblems = initializeSubproblems(obj)
            
            degreesOfFreedom = obj.nShocks-obj.nEqualities-1;
            MaxNumberOfActiveSR = min( obj.nInequalities, degreesOfFreedom);
            % SR stands for Sign Restrictions 
            
            iSubProblem = 1; 
            activeSet = false(obj.nInequalities,1);
            subProblems = subproblemsGivenActiveSet(activeSet);
            
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
 
        
    end
    
    methods (Static)
        function indicesOfActiveSRForSubsets = getAllsubsetsSize_k_from_N(N,k)
            if k==0  % no active restrictions
                indicesOfActiveSRForSubsets = []; 
            else
                indicesOfActiveSRForSubsets = combnk(1: N, k); 
            end
            
        end
    end
end

