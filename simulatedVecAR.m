classdef simulatedVecAR < VecAR
    %simulatedVecAR Summary of this class goes here
    %   Detailed explanation goes here
    

    properties (Access = public )
        simulationSeed = []; % seedMC: seed for a given MC sample
        theta = [];
        thetaCov = [];
        T = [];
    end
    
    methods
        function obj = simulatedVecAR( simulationSeed)
            if nargin~=0
                obj.simulationSeed = simulationSeed;
            end
        end
        function obj = generateThetaFromNormal(obj,design)

            obj.thetaCov  = design.getOmega;
            thetaMean = design.getTheta;
            obj.T = design.getT;
            
            rng('default');
            rng(obj.seedMC,'twister');
            d = size(thetaMean,1);
            stdNormal = randn(d,1);
            
            obj.theta  =  thetaMean + ((obj.thetaCov)^(1/2)/sqrt(obj.T)) * stdNormal;
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

    end
    
end

