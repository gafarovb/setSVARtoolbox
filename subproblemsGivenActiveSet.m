classdef subproblemsGivenActiveSet
    %SUBPROBLEMSFORSVARBOUNDS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        activeSet;
        fullProblem;
    end
    
    methods
        function obj = subproblemsGivenActiveSet(activeSet, objFullProblem)
            obj.activeSet = activeSet;
            obj.fullProblem = fullProbleml;
        end
        
    end
    
end

% -------------------------------------------------------------------------
%% Define nested functions

function [BoundLocal,BoundLocalNeg,Hlocal,valuleLocal] = subproblemBounds(activeSR,bindingSR)
% this subfunction computes the upper bounds BoundU and the corresponding
% argmaxes HU given that Zactive constraints are active and sign constraints Z are imposed


%% initialize local variables
activeZ  = [ZR; SR(activeSR,:)];
inactZ    = SR(~activeSR,:);
BoundLocal    = zeros(nShocks,hori+1);
BoundLocalNeg = zeros(nShocks,hori+1);
valuleLocal    = zeros(nShocks,hori+1);
Hlocal = zeros(nShocks,nShocks,hori+1);


%% 'effective' covariance matrix under active constraints
%  Compare with Lemma 1
if ~isempty(activeZ) % if there are active zero constraints
    M = eye(nShocks) - ( (sigmaSqrt * activeZ')*((activeZ * Sigma * activeZ')^(-1))* (activeZ * sigmaSqrt)) ; % compare with  Proposition 1
    sigmaM = sigmaSqrt * M * sigmaSqrt;   % "Effective" covariance matrix given the set Zactive
else % if there are no active or zero constraints, follow the formula for the  unresrticted case
    sigmaM = Sigma ;
end

for iHorizon=1:hori+1  % loop over IRF horizon
    
    %% bounds under active constraints
    Htmp = sigmaM * reducedIRFaux(:,:,iHorizon)';
    % Compute max value under active constraints
    % - abs is used to avoid complex numbers with Imaginary part equal to Numerical zero
    boundtmp = abs((diag(reducedIRFaux(:,:,iHorizon) * Htmp)).^.5);
    % Compute argmax under active constraints
    %  - bsxfun divides by zeros without warning
    Hlocal(:,:,iHorizon) = bsxfun(@rdivide,Htmp,boundtmp');
    BoundLocal(:,iHorizon)    =   boundtmp ;
    BoundLocalNeg(:,iHorizon) = - boundtmp ;
    valuleLocal(:,iHorizon)   =   boundtmp ;
    
    %% Check whether signs restrictions are satisfied. Compare with Proposition 1
    noSelfPenalty  = true(nShocks,1)  ;
    slackness = zeros(size(inactZ,1),nShocks); % this matrix helps to avoid "selfpenalization"
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
    boundsToPenalize =    (ones(nShocks,1) - ( all(inactZ*Hlocal(:,:,iHorizon)+slackness>-smallNum,1)' & noSelfPenalty ) );
    boundsToPenalizeNeg = (ones(nShocks,1) - ( all(inactZ*Hlocal(:,:,iHorizon)-slackness< smallNum,1)' & noSelfPenalty ) );
    
    
    %% impose penalty for the violation of feasibility constraint
    if ~isempty(inactZ)
        BoundLocal(:,iHorizon)    =   BoundLocal(:,iHorizon) - largeNum * boundsToPenalize;
        BoundLocalNeg(:,iHorizon) =BoundLocalNeg(:,iHorizon) - largeNum * boundsToPenalizeNeg;
    end
    
end


end
% -------------------------------------------------------------------------


