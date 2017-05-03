classdef SVARMomentInequalitiesFramework < handle
    %SVARMomentInequalitiesFramework Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        SVARfacade = [];
        momentInequalities  = [];
        cache;
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
            if isfield(obj.cache, 'cholSigma')
                cholSigma = obj.cache.cholSigma;
            else
                Sigma = obj.SVARfacade.getSigma;
                cholSigma = chol(Sigma);
                obj.cache.cholSigma =  cholSigma;
            end
            xStar = cholSigma * gridpoint;
        end
        function slack = computeInequlaitySlackAtSphericalGridPoint(obj,gridpoint)
            xStar = computeFirstColumnOfSigmaSqrtAtSphericalGridPoint(obj,gridpoint);

            if isfield(obj.cache, 'inequlitiesMat')
                inequlitiesMat = obj.cache.inequlitiesMat;
            else
                linearConstraintsAndDerivatives =  obj.SVARfacade.getLinearConstraintsAndDerivatives;
                inequlitiesMat = linearConstraintsAndDerivatives.SR_nInequalities_sh;
            end
            slack = inequlitiesMat * xStar;
        end
        function residual = computeEqualityResidualAtSphericalGridPoint(obj,gridpoint)
            xStar = computeFirstColumnOfSigmaSqrtAtSphericalGridPoint(obj,gridpoint);
            if isfield(obj.cache, 'equlitiesMat')
                equlitiesMat = obj.cache.equlitiesMat;
            else
                linearConstraintsAndDerivatives =  obj.SVARfacade.getLinearConstraintsAndDerivatives;
                equlitiesMat = linearConstraintsAndDerivatives.ZR_nEqualities_sh;
            end
            residual = equlitiesMat * xStar;
        end
        function valuesIRF = computeObjectiveFunctionsAtSphericalGridPoint(obj,gridpoint)
            xStar = computeFirstColumnOfSigmaSqrtAtSphericalGridPoint(obj,gridpoint);
            if isfield(obj.cache, 'VMA_ts_sh_ho')
                VMA_ts_sh_ho= obj.cache.VMA_ts_sh_ho;
            else
                objectiveFunctionsAndDerivatives = obj.SVARfacade.getIRFObjectiveFunctions;
                VMA_ts_sh_ho = objectiveFunctionsAndDerivatives.VMA_ts_sh_ho;
              
                obj.cache.VMA_ts_sh_ho = VMA_ts_sh_ho;
            end
            valuesIRF = tensorOperations.convWithVector( VMA_ts_sh_ho, 2, xStar);
        end
        function confidenceBounds = twosidedIRFCSbonferroni(obj,level)
            setMomentInequalities(obj);
            
            [confidenceBoundsLow,confidenceBoundsUp] = obj.momentInequalities.computeBonferroniCS(level);
            
            confidenceBoundsLow = obj.SVARfacade.enforceIRFRestrictions(confidenceBoundsLow) ;
            confidenceBoundsUp = obj.SVARfacade.enforceIRFRestrictions(confidenceBoundsUp) ;
             
            confidenceBounds = [confidenceBoundsUp;confidenceBoundsLow];
            
        end
    end
    
end

