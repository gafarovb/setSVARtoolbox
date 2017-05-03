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
        cache;
    end
    
    methods  % constructors
        function obj = SVAR(VecARmodel, ID)
            if nargin > 0
                obj.VecARmodel = VecARmodel;
            else
                obj.VecARmodel = estimatedVecAR; % create a reduced form VAR model
            end
            
            config = obj.getConfig;
            if nargin <2
                warning('estimatedVecAR:identification','Identification assumptions are not provided. Using the default scheme')
                restMat = load(config.assumptionsFilename);
                ID = IDassumptions( restMat);
            end

            obj.label = config.SVARlabel;
            obj.ID = ID;
            
            loadMIframework(obj);
            loadAnalyticFramwork(obj);
        end
        function Samples = generateSamplesFromAsymptoticDistribution(obj,nSimulations)
            rng('default');
            config = obj.getConfig;
            rng(config.masterSeed,'twister');
            seedVector = randi( 1e7, nSimulations); % controls random number generation.
            evalc('Samples(nSimulations) = SVAR;'); % preallocate memory; evalc suppresses the warnings
            
            waitBarWindow = waitBarCustomized(  nSimulations);
            waitBarWindow.setMessage('Generating bootsrap samples...');
            
            for i = 1 : nSimulations
                sampleVecAR = simulatedVecAR(seedVector(i), obj);
                Samples(i) = SVAR(sampleVecAR, obj.ID);
                waitBarWindow.showProgress( i);
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
        function  unitsOfMeasurement = getUnitsOfMeasurement(obj)
            unitsOfMeasurement =  obj.VecARmodel.getUnitsOfMeasurement;
        end 
        function tsDescr = getTSDescription(obj)
           tsDescr =  [obj.getNamesOfTS;obj.getUnitsOfMeasurement];
        end
        
        function shockLabel = getShockLabel(obj)
            shockLabel = obj.ID.shockLabel;
        end
        function MaxHorizons = getMaxHorizons(obj)
            MaxHorizons = obj.VecARmodel.getMaxHorizons;
        end
        function OmegaT = getCovarianceOfThetaT(obj)
            OmegaT = obj.VecARmodel.getCovarianceOfThetaT;
        end
        function restMat = getRestMat(obj)
            restMat = obj.ID.getRestMat;
        end
        function linearConstraintsAndDerivatives = getLinearConstraintsAndDerivatives(obj)
            if isfield(obj.cache,'linearConstraintsAndDerivatives')
                linearConstraintsAndDerivatives = obj.cache.linearConstraintsAndDerivatives;
            else
                linearConstraintsAndDerivatives = obj.ID.getLinearConstraintsAndDerivatives(obj.VecARmodel);
                obj.cache.linearConstraintsAndDerivatives= linearConstraintsAndDerivatives;
            end
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
    
    
    methods (Access = public)
        function IRFCS = IRFtwoSidedCS(obj,level,type)
            lineBreak = char(10);

            if nargin<2
                warning(['Please provide confidence level.' lineBreak ...
                    'Using level = 0.68  '])
                level = 0.68;
            end
            isCorrectLevelRange = (level>0) & (level<1);
            assert( isCorrectLevelRange, 'Please provide confidence level 0<p<1')
     
            
            if nargin<3
                warning(['Please provide type "Analytic" or "MI_Bonferroni".' lineBreak ...
                    ' Running Analytic option '])
                type = 'Analytic';
            end
            
            switch type
                case 'Analytic'
                    disp(['Computing the two-sided CS with p=' num2str(level) ' using the analytic algorithm'  lineBreak ...
                        'from Gafarov, B., M. Meier, and J. L. Montiel Olea (2017) ' lineBreak ...
                        '"DELTA-METHOD INFERENCE FOR A CLASS OF SET-IDENTIFIED SVARS"  ']);
                    IRFCS = obj.analytic.twoSidedIRFCS(level);
                case 'MI_Bonferroni'
                    disp(['Computing the two-sided CS  p=' num2str(level) ' using the Bonferroni-type algorithm'  lineBreak ...
                        'from   Moon, H. R., F. Schorfheide, and E. Granziera (2013) ' lineBreak ...
                        '"Inference for VARs Identified with Sign Restrictions"  ']);
                    IRFCS = obj.MIframework.twosidedIRFCSbonferroni( level);
                otherwise
                warning(['Please provide type "Analytic" or "MI_Bonferroni".' lineBreak ...
                    ' Running Analytic option ']);
                IRFCS = IRFtwoSidedCS(obj,level,'Analytic');
            end
            
        end
        function IRF = estimatedIdentifiedSet(obj)
            IRF = obj.analytic.identifiedSet;
        end
    end
    
end


