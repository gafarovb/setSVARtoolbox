classdef SVAR < handle
    
    %%        SVAR class describes a set-identified SVAR model and offers tools
    %       to construct point estimates and conduct inteference on the impulse
    %       response functions (IRF).
    %
    %       This version: 4.6.2017
    %       Last edited by Bulat Gafarov
    %%
    
    
    properties (Access = public)
        label = 'Unknown'; % Model label, e.g. MSG
        MIframework; 
        analytic;
    end
    
    properties (Access = private)
         VecARmodel  = [];
        ID = [] ; %% An object with restrictions
    end
    
    methods  % constructors
        function obj = SVAR(VecARmodel,restMat)
            if nargin > 0
                obj.VecARmodel = VecARmodel;
            else
                obj.VecARmodel = estimatedVecAR; % create a reduced form VAR model
            end

            config = obj.getConfig;
            if nargin <2
                restMat = load(config.assumptionsFilename);
            end
            
            obj.label = config.SVARlabel;
            obj.ID = IDassumptions( restMat);  
            loadMIframework(obj);
            loadAnalyticFramwork(obj);
        end
        function Samples = generateSamplesFromAsymptoticDistribution(obj,nSimulations)
            rng('default');
            config = obj.getConfig;
            rng(config.masterSeed,'twister');
            seedVector = randi( 1e7, nSimulations); % controls random number generation.
            Samples(nSimulations)=SVAR; % preallocate memory
            for i = 1 : nSimulations
                sampleVecAR = simulatedVecAR(seedVector(i), obj);
                Samples(i) = SVAR(sampleVecAR);
            end
            
        end
        function loadMIframework(obj)
            obj.MIframework = SVARMomentInequalitiesFramework(obj);
        end
        function loadAnalyticFramwork(obj)
            obj.analytic = SVARanalyticFramework(obj);
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
            MaxHorizons = config.nNoncontemoraneousHorizons;
            
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
    end
    

    
    
end


