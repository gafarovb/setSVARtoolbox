classdef simulatedVecAR < VecAR
    %simulatedVecAR Summary of this class goes here
    %   Detailed explanation goes here
    

    properties (Access = public )
        seedMC = []; % seedMC: seed for a given MC sample
        theta = [];
        thetaCov = [];
        T = [];
    end
    
    methods
        function obj = simulatedVecAR( seedMC)
            if nargin~=0
                obj.seedMC = seedMC;
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
    end
    
end

