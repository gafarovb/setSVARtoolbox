classdef SVARanalyticFramework < handle
    %SVARANALYTICFRAMEWORK Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        SVARfacade = [];
        optProblems  ;
    end
    
    methods
        function obj = SVARanalyticFramework(SVARobj)
            obj.SVARfacade = SVARobj;
            obj.optProblems = [];
        end
        function setOptimizationProblems(obj)
            if isempty( obj.optProblems)
                obj.optProblems = optimizationProblems( obj.SVARfacade);
            end
        end
    end
    methods (Access = public)
        function stdIRFcollection = asymptoticStdDeviations(obj)
            stdMat = obj.optProblems.getWorstCaseStdMat;
            stdIRFcollection = IRFcollection(stdMat, obj.SVARfacade.getNamesOfTS, 'standard deviatons') ;
        end
        function maxIRFcollection = onesidedUpperIRFHat(obj)
            setOptimizationProblems(obj);
            maxBoundsMatrix = obj.optProblems.getMaxBounds;
            maxIRFcollection = IRFcollection(maxBoundsMatrix, obj.SVARfacade.getNamesOfTS, 'max point estimates') ;
            maxIRFcollection = obj.SVARfacade.enforceIRFRestrictions(maxIRFcollection) ;
        end
        function minIRFcollection = onesidedLowerIRFHat(obj)
            setOptimizationProblems(obj);
            minBoundsMatrix = obj.optProblems.getMinBounds;
            minIRFcollection = IRFcollection(minBoundsMatrix, obj.SVARfacade.getNamesOfTS, 'min point estimates') ;
            minIRFcollection = obj.SVARfacade.enforceIRFRestrictions(minIRFcollection) ;
        end
        
        function idSet = identifiedSet(obj)
            idSet = [ onesidedUpperIRFHat(obj) ; onesidedLowerIRFHat(obj)];
        end
        function confidenceBounds = onesidedUpperIRFCS(obj,level)
            setOptimizationProblems(obj);
            
            pointEstimates = obj.onesidedUpperIRFHat;
            
            stdIRF =  obj.asymptoticStdDeviations ;
            
            confidenceBounds = pointEstimates +  norminv(level,0,1) * stdIRF  ;
            
            confidenceBounds = obj.SVARfacade.enforceIRFRestrictions( confidenceBounds) ;
            
            confidenceBounds = confidenceBounds.setLabel(['analytic upper one-sided CS with p=',num2str(level)]);
            confidenceBounds = confidenceBounds.setMarker(':');
        end
        function confidenceBounds = onesidedLowerIRFCS(obj,level)
            setOptimizationProblems(obj);
            
            pointEstimates = obj.onesidedLowerIRFHat;
            
            stdIRF =  obj.asymptoticStdDeviations ;
            
            confidenceBounds = pointEstimates - norminv(level,0,1) * stdIRF  ;
            
            confidenceBounds = obj.SVARfacade.enforceIRFRestrictions(confidenceBounds) ;
            
            confidenceBounds = confidenceBounds.setLabel(['analytic lower one-sided CS with p=',num2str(level)]);
            confidenceBounds = confidenceBounds.setMarker(':');
            
        end
        
        function confidenceBounds = twoSidedIRFCS(obj,level)
            levelOneSided = 1- ( 1 - level)/2;
            confidenceBounds = [ obj.onesidedUpperIRFCS(levelOneSided) ; obj.onesidedLowerIRFCS(levelOneSided)];
            confidenceBounds = confidenceBounds.setLabel(['analytic two-sided CS with p=',num2str(level)]);
        end
        
    end
    
    
end

