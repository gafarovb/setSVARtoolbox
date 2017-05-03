classdef SVARanalyticFramework < handle
    %SVARANALYTICFRAMEWORK Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        SVARfacade = [];
        optProblems  ;
    end
    
    methods
        function obj = SVARanalyticFramework( SVARobj)
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
            stdIRFcollection = IRFcollection(stdMat, obj.SVARfacade.getTSDescription, descriptionStd(obj)) ;
        end
        function maxIRFcollection = onesidedUpperIRFHat(obj)
            setOptimizationProblems(obj);
            maxBoundsMatrix = obj.optProblems.getMaxBounds;
            descr = descriptionPointEstimates(obj);
            descr.legend = ['Max' descr.legend ];
            maxIRFcollection = IRFcollection(maxBoundsMatrix, obj.SVARfacade.getTSDescription, descr) ;
            maxIRFcollection = obj.SVARfacade.enforceIRFRestrictions(maxIRFcollection) ;
        end
        function minIRFcollection = onesidedLowerIRFHat(obj)
            setOptimizationProblems(obj);
            minBoundsMatrix = obj.optProblems.getMinBounds;
            descr = descriptionPointEstimates(obj);
            descr.legend = ['Min' descr.legend ];
            minIRFcollection = IRFcollection(minBoundsMatrix, obj.SVARfacade.getTSDescription, descr) ;
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
     
            confidenceBounds = confidenceBounds.setDescriptionField( 'type', 'CS');
            confidenceBounds = confidenceBounds.setDescriptionField( 'legend',[ 'Analytic upper one-sided CS with p=',num2str(level)]);
        end
        function confidenceBounds = onesidedLowerIRFCS(obj,level)
            setOptimizationProblems(obj);
            
            pointEstimates = obj.onesidedLowerIRFHat;
            
            stdIRF =  obj.asymptoticStdDeviations ;
            
            confidenceBounds = pointEstimates - norminv(level,0,1) * stdIRF  ;
            
            confidenceBounds = obj.SVARfacade.enforceIRFRestrictions(confidenceBounds) ;
            
            confidenceBounds = confidenceBounds.setDescriptionField( 'legend',[ 'Analytic lower one-sided CS with p=',num2str(level)]);
            confidenceBounds = confidenceBounds.setDescriptionField( 'type', 'CS');
            
        end
        
        function confidenceBounds = twoSidedIRFCS(obj,level)
            levelOneSided = 1- ( 1 - level)/2;
            confidenceBounds = [ obj.onesidedUpperIRFCS(levelOneSided) ; obj.onesidedLowerIRFCS(levelOneSided)];
            confidenceBounds = confidenceBounds.setDescriptionField( 'legend', ['Analytic two-sided CS with p=',num2str(level)]);
        end
        
        
        
        function  IRFDescription = descriptionPointEstimates(obj)
            SVARobj = obj.SVARfacade;
            
            IRFDescription = struct('tag','noTag',... % can be used in filenames
                'legend','Point estimates',...
                'shock', SVARobj.getShockLabel,...
                'SVARmodel', SVARobj.label,...
                'type' , 'pointEstimates');
            
        end
        
                
        function  IRFDescription = descriptionStd(obj)
                IRFDescription = descriptionPointEstimates(obj);
                IRFDescription.legend = 'Standard Deviations';
                IRFDescription.type = 'std';
                
        end
    end
    
    
end

