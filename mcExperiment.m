classdef mcExperiment
    %MCSVAR Summary of this class goes here
    % This class describes a Monte Carlo exmperiment
    %   Detailed explanation goes here
    
    properties
        design;                  % SVAR object with design
        Samples;                 % Array of SVAR objects with generated samples
    end
    
    properties (Access = private)
        config;
    end
    
    methods
        function obj = mcExperiment(SVARobj)
            obj.design = SVARobj;
            obj.config = obj.design.getConfig;
           
            rng('default');
            rng(obj.config.masterSeed,'twister');
            seedVector = randi(1e7,obj.config.MaxSimulations); % controls random number generation.
            for iMC = 1:obj.config.MaxSimulations
                sampleVecAR = simulatedVecAR(seedVector(iMC), obj.design); 
                Samples(iMC) = SVAR(sampleVecAR);
            end
            obj.Samples = Samples;
        end
        function coverage = testCS(obj)
            % test coverage ??
            coverage = [];
        end
    end
    
end

