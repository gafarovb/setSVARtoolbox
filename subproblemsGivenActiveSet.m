classdef subproblemsGivenActiveSet
    %SUBPROBLEMSFORSVARBOUNDS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        linearEqualityConstraints;
        inactiveInequlities =  [];
        fullProblem;
        activeSet;
        positiveKKTpoints_sh_ts_ho; % indexed : (reduced form ) shock (component) - time series - horizon
        positiveKKTvalues;
        penaltyForPositiveKKT;
        penaltyForNegativeKKT;
        maxBounds;
        minBounds;
    end
    
    methods
        function obj = subproblemsGivenActiveSet(activeSet, objFullProblem)
            obj.activeSet = activeSet;
            
            obj.fullProblem = objFullProblem;
            
            fullSetOfLinearConstraints = objFullProblem.getLinearConstraints;
            
            obj.linearEqualityConstraints  = [fullSetOfLinearConstraints.ZR; ...
                fullSetOfLinearConstraints.SR(activeSet,:)];
            
            obj.inactiveInequlities  = fullSetOfLinearConstraints.SR(~activeSet,:);
            
            obj = computeKKTpointsAndValues(obj);
            obj = computeFeasibilityPenalty(obj);
            obj = computePenalizedMaximum (obj);
            obj = computePenalizedMinimum (obj);
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
        function obj = computeKKTpointsAndValues(obj)
            
            sigmaM = obj.getProjectedSigma;
            objectiveFunctions = obj.fullProblem.getObjectiveFunctions;
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
    end
    
end

