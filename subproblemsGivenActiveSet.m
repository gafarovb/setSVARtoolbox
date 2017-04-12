classdef subproblemsGivenActiveSet < handle  
    %SUBPROBLEMSFORSVARBOUNDS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        linearEqualityConstraints;
        derivativesOfConstraints_nSZ_sh_dAL;
        inactiveInequlities =  [];
        fullProblem;
        activeSet;
        positiveKKTpoints_sh_ts_ho; % indexed : (reduced form ) shock (component) - time series - horizon
        positiveKKTvalues;
        penaltyForPositiveKKT;
        penaltyForNegativeKKT;
    end
    
    properties (Access = public)
        maxBounds;
        minBounds;
    end
    
    
    methods
        function obj = subproblemsGivenActiveSet(activeSet, objFullProblem)
            obj.activeSet = activeSet;  
            obj.fullProblem = objFullProblem;
            
            obj.setActiveConstraintsAndDerivatives;
            obj.computeKKTpointsAndValues ;
            obj.computeFeasibilityPenalty ;
            obj.computePenalizedMaximum ;
            obj.computePenalizedMinimum ;
         end
        
        function obj = setActiveConstraintsAndDerivatives(obj)
            includedInequalities = obj.activeSet;
            
            fullProblemInterface = obj.fullProblem.getLinearConstraintsAndDerivatives;
            
            obj.linearEqualityConstraints  = [fullProblemInterface.ZR_nEqualities_sh; ...
                fullProblemInterface.SR_nInequalities_sh(includedInequalities,:)];
            
            obj.derivativesOfConstraints_nSZ_sh_dAL  = cat( 1, fullProblemInterface.GZR_nEqualities_sh_dAL, ...
                                                 fullProblemInterface.GSR_nInequalities_sh_dAL(includedInequalities,:,:));
            
            obj.inactiveInequlities  = fullProblemInterface.SR_nInequalities_sh(~includedInequalities,:);
        end
        function SigmaProjected = getProjectedSigma(obj)
            %% 'effective' covariance matrix under active constraints
            %  Compare with Lemma 1
            nShocks = size( obj.linearEqualityConstraints , 2);
            sigmaSqrt = obj.fullProblem.getSQRTSigma;
            
            if ~isempty(obj.linearEqualityConstraints) % if there are active zero constraints
                Z = obj.linearEqualityConstraints;
                M = eye(nShocks) - ( (sigmaSqrt * Z')*((Z * sigmaSqrt*sigmaSqrt * Z')^(-1))* (Z * sigmaSqrt)) ; % compare with  Proposition 1
                SigmaProjected = sigmaSqrt * M * sigmaSqrt;   % "Effective" covariance matrix given the set Zactive
            else % if there are no active or zero constraints, follow the formula for the  unresrticted case
                SigmaProjected = sigmaSqrt*sigmaSqrt ;
            end
        end
        function  obj = computeKKTpointsAndValues(obj)
            
            sigmaM = obj.getProjectedSigma;
            objectiveFunctionsAndDerivatives = obj.fullProblem.getObjectiveFunctions;
            objectiveFunctions =  objectiveFunctionsAndDerivatives.VMA_ts_sh_ho;
            [nShocks,~,maxHorizon] = size(objectiveFunctions);
            
            obj.positiveKKTpoints_sh_ts_ho = 0 * objectiveFunctions;
            obj.positiveKKTvalues = zeros(nShocks,maxHorizon);
            
            for aHorizon = 1:maxHorizon
                
                projectedObjectiveFunctionForaHorizon = sigmaM * objectiveFunctions(:,:,aHorizon)';
                
                % Compute max value under active constraints
                % - abs is used to avoid complex numbers with Imaginary part equal to Numerical zero
                obj.positiveKKTvalues(:,aHorizon) = abs((diag(objectiveFunctions(:,:,aHorizon) * projectedObjectiveFunctionForaHorizon)).^.5);
                
                % Compute argmax under active constraints
                %  - bsxfun divides by zeros without warning
                obj.positiveKKTpoints_sh_ts_ho(:,:,aHorizon) = bsxfun(@rdivide, projectedObjectiveFunctionForaHorizon, obj.positiveKKTvalues(:,aHorizon)');
                
            end
        end
        function  obj = computeFeasibilityPenalty(obj)
            
            config = obj.fullProblem.getConfig;
            smallNum = config.smallNumber;
            largeNum = config.largeNumber;
            
            KKT = obj.positiveKKTpoints_sh_ts_ho;
            [nShocks,~,MaxHorizons] = size(KKT);
            feasibileForPositiveKKT = zeros (nShocks,MaxHorizons);
            feasibileForNegativeKKT = zeros (nShocks,MaxHorizons);
            
            nIncative = size(obj.inactiveInequlities ,1);
            slackness = zeros (nIncative,nShocks,MaxHorizons);
            
            for iHorizon = 1:MaxHorizons
                if nIncative>0
                    slackness(:,:,iHorizon) = obj.inactiveInequlities * KKT(:,:,iHorizon);
                end
                feasibileForPositiveKKT(:,iHorizon) = all(   slackness(:,:,iHorizon) > - smallNum,1)';
                feasibileForNegativeKKT(:,iHorizon) = all( - slackness(:,:,iHorizon) > - smallNum,1)';
            end
            obj.penaltyForPositiveKKT =  largeNum * (~feasibileForPositiveKKT) ;
            obj.penaltyForNegativeKKT =  largeNum * (~feasibileForNegativeKKT) ;
        end
        function  obj = computePenalizedMaximum (obj)
            
            maxBoundsPositive =   obj.positiveKKTvalues - obj.penaltyForPositiveKKT;
            maxBoundsNegative = - obj.positiveKKTvalues - obj.penaltyForNegativeKKT;
            
            obj.maxBounds = max( maxBoundsPositive,maxBoundsNegative);
        end
        function  obj = computePenalizedMinimum (obj)
            minBoundsPositive =   obj.positiveKKTvalues + obj.penaltyForPositiveKKT;
            minBoundsNegative = - obj.positiveKKTvalues + obj.penaltyForNegativeKKT;
            
            obj.minBounds = min( minBoundsPositive,minBoundsNegative );
        end
        
        function  stdMat = asymptoticStandardDeviationForActiveSet(obj)

            OmegaT = obj.fullProblem.getCovarianceOfThetaT;
            Sigma = obj.fullProblem.getSigma;
            VMA  =     obj.fullProblem.getObjectiveFunctions;
            C_ts_sh_ho = VMA.VMA_ts_sh_ho;
            G_ts_sh_ho_dAL = VMA.DVMA_ts_sh_ho_dAL;
            [nShocks,~,maxHorizons] = size(C_ts_sh_ho);
            
            xStar_sh_ts_ho = obj.positiveKKTpoints_sh_ts_ho; 
            val_ts_ho = obj.positiveKKTvalues;
            Zactive = obj.linearEqualityConstraints;
            Gactive_nSZ_sh_dAL = obj.derivativesOfConstraints_nSZ_sh_dAL;
            
            dAL = size(Gactive_nSZ_sh_dAL,3);
            restrictionsDerivativeAL = zeros(nShocks,maxHorizons,dAL) ;
            objFunctionDerivativeAL = zeros(nShocks,maxHorizons,dAL) ;
            
            stdMat = zeros(nShocks,maxHorizons);
            
            
            vechFromVec = converter.getVechFromVec(nShocks);   
            for ts = 1:nShocks
                for ho = 1 : maxHorizons
                    wStarForTsHo_nSZ = (Zactive * Sigma *  Zactive')\Zactive * Sigma * C_ts_sh_ho(ts,:,ho)';
                    
                    for iAL = 1:dAL
                        restrictionsDerivativeAL(ts,ho,iAL) = wStarForTsHo_nSZ' *  Gactive_nSZ_sh_dAL(:,:,iAL) * xStar_sh_ts_ho(:,ts,ho);
                        objFunctionDerivativeAL(ts,ho,iAL)  = G_ts_sh_ho_dAL(ts,:,ho,iAL) * xStar_sh_ts_ho(:,ts,ho);
                    end
                    sigmaInvx = Sigma \ xStar_sh_ts_ho(:,ts,ho);
   
                    GradvechSigma =  vechFromVec * (val_ts_ho(ts,ho)/2) * kron(sigmaInvx, sigmaInvx) ;
            
                    GradAL = objFunctionDerivativeAL(ts,ho,:) - restrictionsDerivativeAL(ts,ho,:);
                    
                    Grad =  [squeeze(GradAL);GradvechSigma]';
                     
                    
                    stdMat(ts,ho) = abs((Grad*OmegaT*Grad')^.5);

                end
            end
            
            
            
            
 
            
            
        end
    end
    
end

