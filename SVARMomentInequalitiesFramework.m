classdef SVARMomentInequalitiesFramework < handle
    %SVARMomentInequalitiesFramework Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        SVARfacade = [];
        momentInequalities  = [];
    end
    
    methods
        function obj = SVARMomentInequalitiesFramework(SVARobj)
            obj.SVARfacade = SVARobj;
        end
        function setMomentInequalities(obj)
            if isempty( obj.momentInequalities)
                obj.momentInequalities = stochasticInequalities(obj.SVARfacade) ;
            end
        end
        function xStar = computeFirstColumnOfSigmaSqrtAtSphericalGridPoint(obj,gridpoint)
            Sigma = obj.SVARfacade.getSigma;
            xStar = chol(Sigma) * gridpoint;
        end
        function slack = computeInequlaitySlackAtSphericalGridPoint(obj,gridpoint)
            xStar = computeFirstColumnOfSigmaSqrtAtSphericalGridPoint(obj,gridpoint);
            linearConstraintsAndDerivatives =  obj.SVARfacade.getLinearConstraintsAndDerivatives;
            slack = linearConstraintsAndDerivatives.SR_nInequalities_sh * xStar;
        end
        function residual = computeEqualityResidualAtSphericalGridPoint(obj,gridpoint)
            xStar = computeFirstColumnOfSigmaSqrtAtSphericalGridPoint(obj,gridpoint);
            linearConstraintsAndDerivatives =  obj.SVARfacade.getLinearConstraintsAndDerivatives;
            residual = linearConstraintsAndDerivatives.ZR_nEqualities_sh * xStar;
        end
        function valuesIRF = computeObjectiveFunctionsAtSphericalGridPoint(obj,gridpoint)
            xStar = computeFirstColumnOfSigmaSqrtAtSphericalGridPoint(obj,gridpoint);
            objectiveFunctionsAndDerivatives = obj.SVARfacade.getIRFObjectiveFunctions;
            VMA_ts_sh_ho = objectiveFunctionsAndDerivatives.VMA_ts_sh_ho;
            valuesIRF = tensorOperations.convWithVector( VMA_ts_sh_ho, 2, xStar);
        end
        function confidenceBounds = twosidedIRFCSbonferroni(obj,level)
            setMomentInequalities(obj);
            
            [confidenceBoundsLow,confidenceBoundsUp] = obj.momentInequalities.computeBonferroniCS(level);
            
            confidenceBoundsLow = obj.SVARfacade.enforceIRFRestrictions(confidenceBoundsLow) ;
            confidenceBoundsLow = confidenceBoundsLow.setLabel(['Bonferroni lower two-sided CS with p=',num2str(level)]);
            
            confidenceBoundsUp = obj.SVARfacade.enforceIRFRestrictions(confidenceBoundsUp) ;
            confidenceBoundsUp = confidenceBoundsUp.setLabel(['Bonferroni upper two-sided CS with p=',num2str(level)]);
            
            confidenceBounds = [confidenceBoundsUp,confidenceBoundsLow];
            
        end
    end
    
end

