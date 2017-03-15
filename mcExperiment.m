classdef mcExperiment
    %MCSVAR Summary of this class goes here
    % This class describes a Monte Carlo exmperiment
    %   Detailed explanation goes here
    
    properties
       MaxSimulations = 1000 ;  % number of Monte Carlo simulations
       seedVector;              % controls random number generation.
       Design;                  % SVAR object with design
       array;
    end
    
    properties (Access = private, Constant)
        masterSeed = 123456789; % Seed that generates seeds for every simulation
    end 
    
    methods 
        function obj = mcExperiment(SVARobj)
            obj.Design = SVARobj;
            rng('default');
            rng(obj.masterSeed,'twister');
            obj.seedVector = randi(1e7,obj.MaxSimulations);
            array=SVARobj; % todo make a proper object array
            disp('FIX ME ! mc Experiment')
            for iMC=1:obj.MaxSimulations
                array = SVARobj.resampleTheta(obj.seedVector(iMC));
            end
            obj.array=array;
        end
        function coverage = testCS(obj)
            % test coverage ??
            coverage = [];
        end
    end
    
end

