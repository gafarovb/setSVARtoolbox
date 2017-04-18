classdef SVAR < handle
    
    %%        SVAR class describes a set-identified SVAR model and offers tools
    %       to construct point estimates and conduct inteference on the impulse
    %       response functions (IRF).
    %
    %       This version: 4.6.2017
    %       Last edited by Bulat Gafarov
    %%
    
    
    properties (Access = public)
        label = 'Unknown' ; % Model label, e.g. MSG
    end
    
    properties (Access = private)
        optimiaztionProblems = [];
        VecARmodel;
        momentInequalities;
        ID = [] ; %% An object with restrictions
    end
    
    methods  % constructors
        function obj = SVAR(VecARmodel)
            if nargin > 0
                obj.VecARmodel = VecARmodel;
            else
                obj.VecARmodel = estimatedVecAR; % create a reduced form VAR model
            end
            config = obj.getConfig;
            obj.label = config.label;
            obj.ID = IDassumptions( config.assumptionsFilename); % read ID assumptions from a file
        end
        function Samples = generateSamplesFromAsymptoticDistribution(obj,nSimulations)
            rng('default');
            config = obj.getConfig;
            rng(config.masterSeed,'twister');
            seedVector = randi( 1e7, nSimulations); % controls random number generation.
            Samples(nSimulations)=SVAR;
            for i = 1 : nSimulations
                sampleVecAR = simulatedVecAR(seedVector(i), obj);
                Samples(i) = SVAR(sampleVecAR);
            end
            
        end
        function setMomentInequalities(obj)
            if isempty( obj.momentInequalities)
                obj.momentInequalities = stochasticInequalities(obj) ;
            end
        end
        function setOptimizationProblems(obj)
            if isempty( obj.optimiaztionProblems)
                obj.optimiaztionProblems = optimizationProblems( obj);
            end
        end
    end
    methods
        function config = getConfig(obj)
            config = obj.VecARmodel.getConfig;
        end
        function nShocks = getN(obj)
            nShocks = obj.VecARmodel.getN;
        end
        function nTimePeriods = getT(obj)
            nTimePeriods = obj.VecARmodel.getT;
        end
        function Sigma = getSigma(obj)
            Sigma = obj.VecARmodel.getSigma;
        end
        function theta = getTheta(obj)
            theta  = obj.VecARmodel.getTheta;
        end
        function names = getNamesOfTS(obj)
            names = obj.VecARmodel.getNames;
        end
        function OmegaT = getCovarianceOfThetaT(obj)
            OmegaT = obj.VecARmodel.getCovarianceOfThetaT;
        end
        function restMat = getRestMat(obj)
            restMat = obj.ID.getRestMat;
        end
        function linearConstraintsAndDerivatives = getLinearConstraintsAndDerivatives(obj)
            linearConstraintsAndDerivatives = obj.ID.getLinearConstraintsAndDerivatives(obj.VecARmodel);
        end
        function objectiveFunctions = getIRFObjectiveFunctions(obj)
            objectiveFunMat = obj.VecARmodel.getIRFObjectiveFunctions;
            objectiveFunDerivatives = obj.VecARmodel.getIRFObjectiveFunctionsDerivatives;
            objectiveFunctions = struct(...
                'VMA_ts_sh_ho',objectiveFunMat,...
                'DVMA_ts_sh_ho_dAL',objectiveFunDerivatives);
        end
    end
    
    methods
        function stdIRFcollection = asymptoticStdDeviations(obj)
            stdMat = obj.optimiaztionProblems.getWorstCaseStdMat;
            stdIRFcollection = IRFcollection(stdMat, obj.VecARmodel.getNames, 'standard deviatons') ;
        end
        function IRFcol = enforceIRFRestrictions(SVARobj,IRFcol)
            %% This function enforces sign and zero restrictions specified in SVARobj
            %  on IRF
            %
            % last modified : March 15 2017
            % by Bulat Gafarov
            
            %% read variables from the input objects
            nShocks = SVARobj.getN;
            restMat = SVARobj.getRestMat;
            nRestrictions = size(restMat,1);
            config = SVARobj.getConfig;
            MaxHorizons = config.MaxHorizons;
            
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
            
            IRFmat = IRFcol.matrixForm;
            IRFmat(restrictedAbove) = min(0,IRFmat(restrictedAbove));
            IRFmat(restrictedBelow) = max(0,IRFmat(restrictedBelow));
            IRFcol = IRFcol.setValues(IRFmat);
        end
        function xStar = computeFirstColumnOfSigmaSqrtAtSphericalGridPoint(obj,gridpoint)
            Sigma = obj.getSigma;
            xStar = chol(Sigma) * gridpoint;
        end
        function slack = computeInequlaitySlackAtSphericalGridPoint(obj,gridpoint)
            xStar = computeFirstColumnOfSigmaSqrtAtSphericalGridPoint(obj,gridpoint);
            linearConstraintsAndDerivatives =  obj.getLinearConstraintsAndDerivatives;
            slack = linearConstraintsAndDerivatives.SR_nInequalities_sh * xStar;
        end
        function residual = computeEqualityResidualAtSphericalGridPoint(obj,gridpoint)
            xStar = computeFirstColumnOfSigmaSqrtAtSphericalGridPoint(obj,gridpoint);
            linearConstraintsAndDerivatives =  obj.getLinearConstraintsAndDerivatives;
            residual = linearConstraintsAndDerivatives.ZR_nEqualities_sh * xStar;
        end
        function valuesIRF = computeObjectiveFunctionsAtSphericalGridPoint(obj,gridpoint)
            xStar = computeFirstColumnOfSigmaSqrtAtSphericalGridPoint(obj,gridpoint);
            objectiveFunctionsAndDerivatives = getIRFObjectiveFunctions(obj);
            VMA_ts_sh_ho = objectiveFunctionsAndDerivatives.VMA_ts_sh_ho;
            valuesIRF = tensorOperations.convWithVector( VMA_ts_sh_ho, 2, xStar);
        end
    end
    
    methods (Access = public)
        function maxIRFcollection = onesidedUpperIRFHatAnalytic(obj)
            setOptimizationProblems(obj);
            maxBoundsMatrix = obj.optimiaztionProblems.getMaxBounds;
            maxIRFcollection = IRFcollection(maxBoundsMatrix, obj.VecARmodel.getNames, 'max point estimates') ;
            maxIRFcollection = obj.enforceIRFRestrictions(maxIRFcollection) ;
        end
        function minIRFcollection = onesidedLowerIRFHatAnalytic(obj)
            setOptimizationProblems(obj);
            minBoundsMatrix = obj.optimiaztionProblems.getMinBounds;
            minIRFcollection = IRFcollection(minBoundsMatrix, obj.VecARmodel.getNames, 'min point estimates') ;
            minIRFcollection = obj.enforceIRFRestrictions(minIRFcollection) ;
        end
        function confidenceBounds = onesidedUpperIRFCSAnalytic(obj,level)
            setOptimizationProblems(obj);
             
            pointEstimates = obj.onesidedUpperIRFHatAnalytic;
            
            stdIRF =  obj.asymptoticStdDeviations ;
            
            confidenceBounds = pointEstimates +  norminv(level,0,1) * stdIRF  ;
            
            confidenceBounds = obj.enforceIRFRestrictions(confidenceBounds) ;
            
            confidenceBounds = confidenceBounds.setLabel(['analytic upper one-sided CS with p=',num2str(level)]);
        end
        function confidenceBounds = onesidedLowerIRFCSAnalytic(obj,level)
            setOptimizationProblems(obj);
                        
            pointEstimates = obj.onesidedLowerIRFHatAnalytic;
            
            stdIRF =  obj.asymptoticStdDeviations ;
            
            confidenceBounds = pointEstimates - norminv(level,0,1) * stdIRF  ;
            
            confidenceBounds = obj.enforceIRFRestrictions(confidenceBounds) ;
            
            confidenceBounds = confidenceBounds.setLabel(['analytic lower one-sided CS with p=',num2str(level)]);
        end
        function confidenceBounds = twosidedIRFCSbonferroni(obj,level)
            setMomentInequalities(obj);
            
            [confidenceBoundsLow,confidenceBoundsUp] = obj.momentInequalities.computeBonferroniCS(level);
            
            confidenceBoundsLow = obj.enforceIRFRestrictions(confidenceBoundsLow) ;
            confidenceBoundsLow = confidenceBoundsLow.setLabel(['Bonferroni lower two-sided CS with p=',num2str(level)]);
            
            confidenceBoundsUp = obj.enforceIRFRestrictions(confidenceBoundsUp) ;
            confidenceBoundsUp = confidenceBoundsUp.setLabel(['Bonferroni upper two-sided CS with p=',num2str(level)]);
            
            confidenceBounds = [confidenceBoundsUp,confidenceBoundsLow];
            
        end
    end
    
    methods 
    end
    
    
end


