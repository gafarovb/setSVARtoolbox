classdef SVAR  
    
    %%        SVAR class describes a set-identified SVAR model and offers tools
    %       to construct point estimates and conduct inteference on the impulse
    %       response functions (IRF).
    %
    %       This version: 4.6.2017
    %       Last edited by Bulat Gafarov
    %%

    
    properties (Access = public)
        label = 'Unknown' ; % Model label, e.g. MSG 
        ID   =[] ; %% An object with restrictions
                momentInequalities;

     end
    
    properties (Access = private)
        optimiaztionProblems;
        VecARmodel;
    end
        
    methods
        function obj = SVAR(VecARmodel)

            if nargin > 0
                obj.VecARmodel = VecARmodel;
            else
                obj.VecARmodel = estimatedVecAR; % create a reduced form VAR model
            end
            config = obj.getConfig;
            obj.label = config.label;
            obj.ID = IDassumptions( config.assumptionsFilename); % read ID assumptions from a file
            
            obj.optimiaztionProblems = optimizationProblems( obj);
        end     
        function config = getConfig(obj)
            config = obj.VecARmodel.getConfig;
        end
        function nShocks = getN(obj)
            nShocks = obj.VecARmodel.getN;
        end
        function nTimePeriods = getT(obj)
            nTimePeriods = obj.VecARmodel.getT;
        end
        function linearConstraintsAndDerivatives = getLinearConstraintsAndDerivatives(obj)
            linearConstraintsAndDerivatives = obj.ID.getLinearConstraintsAndDerivatives(obj.VecARmodel);
        end
        function OmegaT = getCovarianceOfThetaT(obj)
           OmegaT = obj.VecARmodel.getCovarianceOfThetaT; 
        end
        function Sigma = getSigma(obj)
            Sigma = obj.VecARmodel.getSigma;
        end
        function restMat = getRestMat(obj)
            restMat = obj.ID.getRestMat;
        end
        function objectiveFunctions = getIRFObjectiveFunctions(obj)
            objectiveFunMat = obj.VecARmodel.getIRFObjectiveFunctions;
            objectiveFunDerivatives = obj.VecARmodel.getIRFObjectiveFunctionsDerivatives;
            objectiveFunctions = struct(...
            'VMA_ts_sh_ho',objectiveFunMat,...
            'DVMA_ts_sh_ho_dAL',objectiveFunDerivatives);
        end
        function stdIRFcollection = asymptoticStdDeviations(obj)
           stdMat = obj.optimiaztionProblems.getWorstCaseStdMat;
           stdIRFcollection = IRFcollection(stdMat, obj.VecARmodel.getNames, 'standard deviatons') ;
        end
        function IRFcollection = enforceIRFRestrictions(SVARobj,IRFcollection)
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
            
            IRFmat = IRFcollection.matrixForm;
            IRFmat(restrictedAbove) = min(0,IRFmat(restrictedAbove));
            IRFmat(restrictedBelow) = max(0,IRFmat(restrictedBelow));
            IRFcollection = IRFcollection.setValues(IRFmat); 
        end
        function theta = getTheta(obj)
            theta  = obj.VecARmodel.getTheta;
        end
        function names = getNamesOfTS(obj)
            names = obj.VecARmodel.getNames;
        end
        function slack = computeInequlaitySlackAtSphericalGridPoint(obj,gridpoint)
            Sigma = obj.getSigma;
            xStar = chol(Sigma) * gridpoint;
            linearConstraintsAndDerivatives =  obj.getLinearConstraintsAndDerivatives;
            slack = linearConstraintsAndDerivatives.SR_nInequalities_sh * xStar;
        end
        function residual = computeEqualityResidualAtSphericalGridPoint(obj,gridpoint)
            Sigma = obj.getSigma;
            xStar = chol(Sigma) * gridpoint;
            linearConstraintsAndDerivatives =  obj.getLinearConstraintsAndDerivatives;
            residual = linearConstraintsAndDerivatives.ZR_nEqualities_sh * xStar;
        end
        function Samples = generateSamplesFromAsymptoticDistribution(obj,nSimulations)
            rng('default');
            config = obj.getConfig;
            rng(config.masterSeed,'twister');
            seedVector = randi( 1e7, nSimulations); % controls random number generation.
            for i = 1 : nSimulations
                sampleVecAR = simulatedVecAR(seedVector(i), obj);
                Samples(i) = SVAR(sampleVecAR);
            end
             
        end
    end
    
    methods (Access = public)
        function maxIRFcollection = onesidedUpperIRFHatAnalytic(obj)
            maxBoundsMatrix = obj.optimiaztionProblems.getMaxBounds;
            maxIRFcollection = IRFcollection(maxBoundsMatrix, obj.VecARmodel.getNames, 'max point estimates') ;
            maxIRFcollection = obj.enforceIRFRestrictions(maxIRFcollection) ;

        end
        function minIRFcollection = onesidedLowerIRFHatAnalytic(obj)
            minBoundsMatrix = obj.optimiaztionProblems.getMinBounds;
            minIRFcollection = IRFcollection(minBoundsMatrix, obj.VecARmodel.getNames, 'min point estimates') ;
            minIRFcollection = obj.enforceIRFRestrictions(minIRFcollection) ;
        end
        function confidenceBounds = onesidedUpperIRFCSAnalytic(obj,level)
             
            pointEstimates = obj.onesidedUpperIRFHatAnalytic;
            
            stdIRF =  obj.asymptoticStdDeviations ;
            
            confidenceBounds = pointEstimates +  norminv(level,0,1) * stdIRF  ;
            
            confidenceBounds = obj.enforceIRFRestrictions(confidenceBounds) ;
            
            confidenceBounds = confidenceBounds.setLabel(['analytic upper one-sided CS with p=',num2str(level)]);
        end
        function confidenceBounds = onesidedLowerIRFCSAnalytic(obj,level)
            
            pointEstimates = obj.onesidedLowerIRFHatAnalytic;
            
            stdIRF =  obj.asymptoticStdDeviations ;
            
            confidenceBounds = pointEstimates - norminv(level,0,1) * stdIRF  ;
            
            confidenceBounds = obj.enforceIRFRestrictions(confidenceBounds) ;
            
            confidenceBounds = confidenceBounds.setLabel(['analytic lower one-sided CS with p=',num2str(level)]);
        end
        
        
        function obj = setMomentInequalities(obj)
            obj.momentInequalities = stochasticInequalities(obj) ;
        end
    end
    
 
end
 

