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
        function obj = IDassumptions(assumptionsFilename)
            % read the restricitons matrix from a file  'restMat.dat'
            obj.assumptionsMatrixInput = load(assumptionsFilename);
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
        
        function linearConstraintsAndDerivatives = generateLinearConstraints( obj,objVecAR)
            %% -------------------------------------------------------------------------
            % This function generates equality and inequality constraints
            %
            % This version: March 24th, 2017
            % Last edited :  Bulat Gafarov
            % -------------------------------------------------------------------------
            
            nSignRestrictions = obj.countSignRestrictions;
            nZeroRestrictions = obj.countZeroRestrictions;
           
            n =  objVecAR.getN;
            G =  objVecAR.getVMADerivatives;
            VMA    = objVecAR.getVMA_ts_sh_ho;
            Gcum = cumsum(G,3);
            VMAcum = cumsum(VMA,3);
            
            [nG,mG,~] = size( G);
            
            %% allocate memory
            SR = zeros( nSignRestrictions, n); % matrix with sign restrictions
            ZR = zeros( nZeroRestrictions, n); % matrix with zero restrictions
            GSR = zeros( nG, mG, nSignRestrictions); % Derivatives of sign restricted IRF
            GZR = zeros( nG, mG, nZeroRestrictions); % Derivatives of zero restricted IRF
            
            iSR = 0;
            iZR = 0;
            totalNumberOfRestriction = nSignRestrictions+nZeroRestrictions;

            for aRestriction = 1 :totalNumberOfRestriction % loop through all restricions
                
                horizon = obj.getHorizonOfaRestriction( aRestriction);
                isCumulative = obj.aRestrictionIsCumulative( aRestriction);
                tsOfaRestriction = obj.getTsOfaRestriction( aRestriction);
                
                currentIRF = (1-isCumulative) *     VMA( tsOfaRestriction,:, horizon) +...
                    isCumulative  *  VMAcum(tsOfaRestriction,:,horizon);
                
                currentRestrictionDerivative = (1-isCumulative) *    G(:,:, horizon) +...
                    isCumulative  * Gcum(:,:,horizon);
                 
                if  obj.aRestrictionIsEquality( aRestriction)   
                    iZR=iZR+1;
                    ZR(iZR,:)    = currentIRF;
                    GZR(:,:,iZR) = currentRestrictionDerivative;
                else     
                    iSR=iSR+1;
                    SR(iSR,:) = obj.convertArestrictionToGeq(aRestriction) * currentIRF;
                    GSR(:,:,iSR) = obj.convertArestrictionToGeq(aRestriction) * currentRestrictionDerivative;                
                end
            end
            
            linearConstraintsAndDerivatives = struct(...
                'GSR',GSR, ...    % G matrix associated with inequality constraints
                'GZR',GZR, ...    % G matrix associated with   equality constraints
                'SR',SR, ...   % sign restrictions
                'ZR',ZR);      % zero restrictions
            
        end
    end
    
end

