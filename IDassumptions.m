classdef IDassumptions < handle
    
    
    %%      restrictions class describes restrictions on IRF in SVARs
    %
    %       This version: 3.21.2017
    %       Last edited by Bulat Gafarov
    %%
    
    
    properties (Access = private)
        assumptionsMatrixInput;  % Matrix with columns: Var Hor S/Z Cum Shk
    end

    
    methods
        function obj = IDassumptions(restMat)
             obj.assumptionsMatrixInput =restMat ;
        end
        function matrix = getRestMat(obj)
            % read the restricitons matrix
            matrix = obj.assumptionsMatrixInput;
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
        
        
        function horizon = getHorizonOfaRestriction( obj, aRestriction)
            horizon = obj.assumptionsMatrixInput( aRestriction,2) + 1;
        end
        function tsOfaRestriction = getTsOfaRestriction( obj, aRestriction)
            tsOfaRestriction = obj.assumptionsMatrixInput( aRestriction,1);
        end
        function isCumulative = aRestrictionIsCumulative( obj, aRestriction)
            isCumulative = obj.assumptionsMatrixInput( aRestriction,2) + 1;
        end
        function isZeroRestirction = aRestrictionIsEquality( obj, aRestriction)
            isZeroRestirction = (obj.assumptionsMatrixInput( aRestriction,3)==0) ;
        end  
        function signedOne = convertArestrictionToGeq( obj, aRestriction)
            signedOne = obj.assumptionsMatrixInput(aRestriction,3);
        end
        
        function linearConstraintsAndDerivatives = getLinearConstraintsAndDerivatives( obj,objVecAR)
            %% -------------------------------------------------------------------------
            % This function generates equality and inequality constraints
            %
            % This version: March 24th, 2017
            % Last edited :  Bulat Gafarov
            % -------------------------------------------------------------------------
            
            nSignRestrictions = obj.countSignRestrictions;
            nZeroRestrictions = obj.countZeroRestrictions;
           
            nShock =  objVecAR.getN;
            G_ts_sh_ho_dAL =  objVecAR.getVMADerivatives_ts_sh_ho_dAL;
            VMA_ts_sh_ho    = objVecAR.getVMA_ts_sh_ho;
            Gcum_ts_sh_ho_dAL = cumsum(G_ts_sh_ho_dAL,3);
            VMAcum_ts_sh_ho = cumsum(VMA_ts_sh_ho,3);
            dAL = size( G_ts_sh_ho_dAL,4);
            
            %% allocate memory
            SR_nInequalities_sh = zeros( nSignRestrictions, nShock); % matrix with sign restrictions
            ZR_nEqualities_sh = zeros( nZeroRestrictions, nShock); % matrix with zero restrictions
            GSR_sh_nInequalities_dAL = zeros( nShock, nSignRestrictions, dAL); % Derivatives of sign restricted IRF
            GZR_sh_nEqualities_dAL = zeros( nShock, nZeroRestrictions, dAL); % Derivatives of zero restricted IRF
            
            iSR = 0;
            iZR = 0;
            totalNumberOfRestriction = nSignRestrictions+nZeroRestrictions;

            for aRestriction = 1 : totalNumberOfRestriction % loop through all restricions
                
                horizon = obj.getHorizonOfaRestriction( aRestriction);
                isCumulative = obj.aRestrictionIsCumulative( aRestriction);
                tsOfaRestriction = obj.getTsOfaRestriction( aRestriction);
                
                currentIRF = (1-isCumulative) *     VMA_ts_sh_ho( tsOfaRestriction,:, horizon) +...
                    isCumulative  *  VMAcum_ts_sh_ho(tsOfaRestriction,:,horizon);
                
                currentRestrictionDerivative = (1-isCumulative) *    G_ts_sh_ho_dAL(tsOfaRestriction,:, horizon,:) +...
                    isCumulative  * Gcum_ts_sh_ho_dAL(tsOfaRestriction,:,horizon,:);
                 
                if  obj.aRestrictionIsEquality( aRestriction)   
                    iZR=iZR+1;
                    ZR_nEqualities_sh(iZR,:)    = currentIRF;
                    GZR_sh_nEqualities_dAL(:,iZR,:) = currentRestrictionDerivative;
                else     
                    iSR=iSR+1;
                    SR_nInequalities_sh(iSR,:) = obj.convertArestrictionToGeq(aRestriction) * currentIRF;
                    GSR_sh_nInequalities_dAL(:,iSR,:) = obj.convertArestrictionToGeq(aRestriction) * currentRestrictionDerivative;                
                end
            end
 
            GSR_nInequalities_sh_dAL = permute(GSR_sh_nInequalities_dAL,[2,1,3]);
            GZR_nEqualities_sh_dAL = permute(GZR_sh_nEqualities_dAL,[2,1,3]);
            
            linearConstraintsAndDerivatives = struct(...
                'GSR_nInequalities_sh_dAL',GSR_nInequalities_sh_dAL, ...    % G array associated with inequality constraints
                'GZR_nEqualities_sh_dAL',GZR_nEqualities_sh_dAL, ...    % G array associated with   equality constraints
                'SR_nInequalities_sh',SR_nInequalities_sh, ...   % sign restrictions
                'ZR_nEqualities_sh',ZR_nEqualities_sh);      % zero restrictions
            
        end
    end
    
end

