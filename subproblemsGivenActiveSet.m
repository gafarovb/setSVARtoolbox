classdef subproblemsGivenActiveSet
    %SUBPROBLEMSFORSVARBOUNDS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        linearEqualityConstraints;
        inactiveInequlities;
        fullProblem;
        activeSet;
        KKTpoint_sh_ts_ho; % indiced : (reduced form ) shock (component) - time series - horizon
    end
    
    methods
        function obj = subproblemsGivenActiveSet(activeSet, objFullProblem)
            obj.activeSet = activeSet;
            
            obj.fullProblem = objFullProblem;
            
            fullSetOfLinearConstraints = objFullProblem.getLinearConstraints;
            
            obj.linearEqualityConstraints  = [fullSetOfLinearConstraints.ZR; ...
                fullSetOfLinearConstraints.SR(activeSet,:)];
            
            obj.inactiveInequlities    = fullSetOfLinearConstraints.SR(~activeSet,:);
            
            obj.KKTpoint_sh_ts_ho = subproblemKKTpoints(obj);
             subproblemPenalizedMaximum (obj)
            
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
        function KKTpoints = subproblemKKTpoints(obj)
            
            sigmaM = obj.getProjectedSigma;
            objectiveFunctions = obj.fullProblem.getObjectiveFunctions;
            maxHorizon = size(objectiveFunctions,3);
            KKTpoints = 0 * objectiveFunctions;
            
            for aHorizon = 1:maxHorizon
                
                KKTforaHorizon = sigmaM * objectiveFunctions(:,:,aHorizon)';
                
                % Compute max value under active constraints
                % - abs is used to avoid complex numbers with Imaginary part equal to Numerical zero
                KKTvalues = abs((diag(objectiveFunctions(:,:,aHorizon) * KKTforaHorizon)).^.5);
                
                % Compute argmax under active constraints
                %  - bsxfun divides by zeros without warning
                KKTpoints(:,:,aHorizon) = bsxfun(@rdivide,KKTforaHorizon,KKTvalues');
                
            end
        end
        function subPenalizedMax= subproblemPenalizedMaximum (obj)
            
            
            Hlocal = obj.KKTpoint_sh_ts_ho;
            [nShocks,~,maxHorizon] = size(Hlocal);
            
            
            for iHorizon = 1 : maxHorizon  % loop over IRF horizon
                %% Check whether signs restrictions are satisfied. Compare with Proposition 1
                noSelfPenalty  = true(nShocks,1)  ;
                
                nIactiveConstraints  = size( obj.inactiveInequlities,1);
                
                slackness = zeros( nIactiveConstraints, nShocks); % this matrix helps to avoid "selfpenalization"
                
                for iBindingSR= 1 : size(bindingSR,1)
                    
                    if   ( bindingSR(iBindingSR,2)+1 ) == iHorizon
                   
                        tsIndex = bindingSR( iBindingSR,1);                       

                        indexOftheSignRestriction = bindingSR(iBindingSR,4);

                        if activeSR( indexOftheSignRestriction )
                            
                            noSelfPenalty(tsIndex) = false; % this matrix helps to avoid "selfpenalization"
                        else
                            
                            columnWithOne = zeros( size(SR,1),1);
                            
                            columnWithOne(indexOftheSignRestriction) = abs(largeNum);
                            
                            columnWithOne = columnWithOne( ~activeSR,:);
                            
                            slackness( indexOftheSignRestriction , tsIndex )= columnWithOne;
                        end
                        
                    end
                end
                
                % if there are no sign restrictions, the idicator functions are equal to 1
                % iff  some sign restriction is violated by a margin of a numerical error of 1.e-6
                % the corresponding indicator function is equal to 0
                
                boundsToPenalize =    (ones(nShocks,1) - ( all( obj.inactiveInequlities * Hlocal(:,:,iHorizon)+slackness>-smallNum,1)' & noSelfPenalty ) );

                
                
                
                
                
                
                
                
                
                
                
                boundsToPenalizeNeg = (ones(nShocks,1) - ( all( obj.inactiveInequlities * Hlocal(:,:,iHorizon)-slackness< smallNum,1)' & noSelfPenalty ) );
                subPenalizedMax = 0;
            end
        end
    end
    
end
function [BoundLocal,BoundLocalNeg,Hlocal,valuleLocal] = legacys(obj)

%% initialize local variables

%             BoundLocal    = zeros(nShocks,hori+1);
%             BoundLocalNeg = zeros(nShocks,hori+1);
%             valuleLocal    = zeros(nShocks,hori+1);
%             Hlocal = zeros(nShocks,nShocks,hori+1);
%
sigmaM = getProjectedSigma(obj);
objectiveFunctions = obj.fullProblem.getObjectiveFunctions;

for iHorizon=1:hori+1  % loop over IRF horizon
    
    %% bounds under active constraints
    Htmp = sigmaM * objectiveFunctions(:,:,iHorizon)';
    
    % Compute max value under active constraints
    % - abs is used to avoid complex numbers with Imaginary part equal to Numerical zero
    boundtmp = abs((diag(objectiveFunctions(:,:,iHorizon) * Htmp)).^.5);
    
    % Compute argmax under active constraints
    %  - bsxfun divides by zeros without warning
    Hlocal(:,:,iHorizon) = bsxfun(@rdivide,Htmp,boundtmp');
    
    
    BoundLocal(:,iHorizon)    =   boundtmp ;
    BoundLocalNeg(:,iHorizon) = - boundtmp ;
    valuleLocal(:,iHorizon)   =   boundtmp ;
    
    %% Check whether signs restrictions are satisfied. Compare with Proposition 1
    noSelfPenalty  = true(nShocks,1)  ;
    slackness = zeros(size(inactiveInequlities,1),nShocks); % this matrix helps to avoid "selfpenalization"
    for iBindingSR=1:size(bindingSR,1)
        if   (bindingSR(iBindingSR,2)+1)==iHorizon
            if activeSR( bindingSR(iBindingSR,4) )
                noSelfPenalty(bindingSR(iBindingSR,1)) = false; % this matrix helps to avoid "selfpenalization"
            else
                columnWithOne = zeros(size(SR,1),1);
                columnWithOne(bindingSR(iBindingSR,4))=abs(largeNum);
                columnWithOne = columnWithOne(~activeSR,:);
                slackness(:,bindingSR(iBindingSR,1))= columnWithOne;
            end
        end
    end
    
    % if there are no sign restrictions, the idicator functions are equal to 1
    % iff  some sign restriction is violated by a margin of a numerical error of 1.e-6
    % the corresponding indicator function is equal to 0
    boundsToPenalize =    (ones(nShocks,1) - ( all(inactiveInequlities*Hlocal(:,:,iHorizon)+slackness>-smallNum,1)' & noSelfPenalty ) );
    boundsToPenalizeNeg = (ones(nShocks,1) - ( all(inactiveInequlities*Hlocal(:,:,iHorizon)-slackness< smallNum,1)' & noSelfPenalty ) );
    
    
    %% impose penalty for the violation of feasibility constraint
    if ~isempty(inactiveInequlities)
        BoundLocal(:,iHorizon)    =   BoundLocal(:,iHorizon) - largeNum * boundsToPenalize;
        BoundLocalNeg(:,iHorizon) =BoundLocalNeg(:,iHorizon) - largeNum * boundsToPenalizeNeg;
    end
    
end


end
