classdef IDrestrictions<handle
    
    
    %%      restrictions class describes restrictions on IRF in SVARs
    %
    %       This version: 3.21.2017
    %       Last edited by Bulat Gafarov
    %%
    
    
    properties (Access = private)
        restrictionsMatrixInput;  % Matrix with columns: Var Hor S/Z Cum Shk
    end
    
    properties (SetAccess = public)
        selectorSR;    % 'e' vectors of inequality restrictions
        selectorZR;    % 'e' vectors of equality restrictions
        maskColMax;    % IRF to skip
        maskColMin;    % IRF to skip
        negativeSR;    % matrix with indexed negative sign restrictions : Var Horizon (-1) iSignRestriction 
        positiveSR;    % matrix with indexed positive sign restrictions : Var Horizon ( 1) iSignRestriction
    end
    
    methods
        function obj = IDrestrictions(restricitonsFilename)
            % read the restricitons matrix from a file  'restMat.dat'
            obj.restrictionsMatrixInput = load(restricitonsFilename);
        end
        function obj = constructSelectors(obj,nShocks)
            % this function constructs selector matrices corresponding to the restrictions 
            
            %% interface
            nSignRestrictions = countSignRestrictions(obj);
            nZeroRestrictions = countZeroRestrictions(obj);
            MaxHorizons = configFile.MaxHorizons;         
            restMat = obj.getRestMat;
            
            %% allocate memory
            orderOfSR  = zeros(nSignRestrictions,4);  
            obj.selectorSR = zeros(nShocks,nSignRestrictions);
            obj.selectorZR = zeros(nShocks,nZeroRestrictions);
            Emat = eye(nShocks);
            iSR = 1;
            iZR = 1;
            lf  = 0;
            
            for i=1:(nSignRestrictions+nZeroRestrictions) % loop through all restricions
                if  restMat(i,3)==0                      % check if it is a zero restrictions
                    obj.selectorZR(:,iZR) = Emat(:,restMat(i,1));
                    iZR=iZR+1;
                else                                     % assume sign constraints
                    lf=lf+1;
                    orderOfSR(lf,1:3)=restMat(i,1:3);
                    orderOfSR(lf,4)= iSR;
                    obj.selectorSR(:,iSR) = Emat(:,restMat(i,1));
                    iSR=iSR+1;
                end
            end
            
            obj.negativeSR = orderOfSR(orderOfSR(:,3)<0,:);
            obj.positiveSR = orderOfSR(orderOfSR(:,3)>0,:); % matrix with indexed positive sign restrictions
           
            if ~isempty(obj.negativeSR)
                obj.maskColMax  = logical( full(sparse(obj.negativeSR(:,1),        1+ obj.negativeSR(:,2),1,nShocks,MaxHorizons+1)));
            else
                obj.maskColMax  = false(nShocks,MaxHorizons+1)  ;
            end
            
            if ~isempty(obj.positiveSR)
                obj.maskColMin  = logical( full(sparse(obj.positiveSR(:,1),        1+ obj.positiveSR(:,2),1,nShocks,MaxHorizons+1)));
            else
                obj.maskColMin  = false(nShocks,MaxHorizons+1)  ;
            end

            
        end
        function matrix = getRestMat(obj)
            % read the restricitons matrix
            matrix = obj.restrictionsMatrixInput;
        end
        function nSignRestrictions = countSignRestrictions(obj)
            restMat = obj.getRestMat ;
            nSignRestrictions = sum(abs(restMat(:,3))); % number of sign restrictions
        end
        function nZeroRestrictions = countZeroRestrictions(obj)
            restMat = obj.getRestMat;
            nZeroRestrictions = sum(restMat(:,3)==0);   % number of zero restrictions
        end
        function nAcitiveSets      = countActiveSets(obj,n)
            %% compute number of possible combinations of active sign restrictions
            nAcitiveSets = 1;
            nSignRestrictions = obj.countSignRestrictions;
            for iBin = 1:min(n-1,nSignRestrictions)
                nAcitiveSets = nAcitiveSets + nchoosek(nSignRestrictions,iBin);
            end
        end
    end
    
end

