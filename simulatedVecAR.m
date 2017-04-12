classdef simulatedVecAR < VecAR
    %simulatedVecAR Summary of this class goes here
    %   Detailed explanation goes here
    

    properties (Access = public )
        simulationSeed = []; % seedMC: seed for a given MC sample
        theta = [];
        thetaCovT = [];
        nShocks =[];
        namesOfTS=[];
     end
    
    methods
        function obj = simulatedVecAR( simulationSeed , design)
            if nargin~=0
                obj.simulationSeed = simulationSeed;
                obj.namesOfTS = design.getNamesOfTS;
                obj.config = design.getConfig;
                obj.nLags = obj.config.nLags ;
                obj.nShocks = design.getN;
                
                obj = generateThetaFromNormal(obj,design);
                
            end
        end
        function nShocks = getN(obj)
            nShocks  = obj.nShocks;
        end
        function obj = generateThetaFromNormal(obj,design)

            obj.thetaCovT  = design.getCovarianceOfThetaT;
            thetaMean = design.getTheta;
             
            rng('default');
            rng(obj.simulationSeed,'twister');
            d = size(thetaMean,1);
            stdNormal = randn(d,1);
            
            obj.theta  =  thetaMean + ((obj.thetaCovT)^(1/2)) * stdNormal;
        end
        function simVAR   = simulateVAR(AL_n_x_np,T,nMC,etahat)
            %% TODO: Legacy code
            
            % -------------------------------------------------------------------------
            % This function simulate nMC samples of length T based on VAR model with lag polynomial AL
            % and covariance matrix Sigma
            %
            % Inputs:
            % - AL: VAR model coefficients
            % - T: timer series length
            % - nMC: number of bootstrap samples
            % - eta: VAR model residuals
            % Outputs:
            % - simVAR: (T x n x nMC) 3d array of simulated data
            %
            % This version: February 24, 2015
            % -------------------------------------------------------------------------
            
            % To do list
            % - use MA representation to vectorize the simulation for better speed
            
            [n,np] = size(AL_n_x_np);
            lags = np/n;
            burnIn = 1000;
            etahat =  etahat -mean(etahat,2)* ones( 1, T-lags) ;
            
            eta = etahat( randi([1,(T-lags)],nMC*T+burnIn,n) );
            TSLtemp =  zeros( nMC * T + burnIn, n);
            for iT = ( lags+1):( nMC * T + burnIn)
                TSLtemp(iT,:) =  ( AL_n_x_np * reshape( (TSLtemp((iT-1):-1:(iT- lags),:))',[ n *  lags,1]) )' + eta(iT,:);
            end
            
            simVAR = permute ( reshape(TSLtemp((burnIn+1):( nMC* T+burnIn),:),[ T ,  nMC ,  n ]),   [1 3 2]);
            
        end
        function Sigma = getSigma(obj)
            [~,Sigma] =  obj.ALSigmaFromThetaNandP( obj.theta, obj.getN, obj.nLags);
        end
        function AL_n_x_np = getAL_n_x_np(obj)
            [AL_n_x_np,~] =  obj.ALSigmaFromThetaNandP( obj.theta, obj.getN, obj.nLags);
        end
        function thetaCovT = getCovarianceOfThetaT(obj)
            thetaCovT = obj.thetaCovT;
        end
        function theta = getTheta(obj)
            theta = obj.theta;
        end
        function VMA_ts_sh_ho = getVMA_ts_sh_ho(obj)
            AL_n_x_np = obj.getAL_n_x_np;
            hori = obj.config.MaxHorizons;
            VMA_ts_sh_ho = obj.getVMAfromAL( AL_n_x_np, hori);
        end
        function G = getVMADerivatives_ts_sh_ho_dAL(obj)          
            AL_n_x_np = obj.getAL_n_x_np;
            hori = obj.config.MaxHorizons;
            G = obj.getVMAderivatives( AL_n_x_np, hori);
        end
        function names =  getNames(obj)
            names = obj.namesOfTS;
        end
    end
    

end

