classdef configSVAR<handle
    %configSVAR contains specification for SVAR objects
    %   Detailed explanation goes here
    
    properties (Access = public,Constant )
                %% Specify shock
        % solve for IRF bounds of specified shock 
        % (by default, compute bounds for all shocks under gridsearch algorithm)
    masterSeed = 123456789; % Seed that generates seeds for every simulation
    MaxSimulations = 1000 ;  % number of Monte Carlo simulations
    
    shock = 1; 
%% further model options
    coverage = 0; % 0=off, 1=on: compute MC coverage frequency
    cluster = 0;  % 0=off, 1=on: use cluster to compute MC coverage frequency
    caliProj = 0; % 0=off, 1=on: 'fake' calibration of Projection CS, 2=on: calibration using IRF rotations 
    sobol = 0;    % 0=off, 1=on: use Sobol sequences instead of random numbers
    bootstrap = 0; % 0=draw reduced-form parameters, 1=bootstrap data and re-estimate reduced-form parameters (instead of drawing them directly)
    dif_calibration = 0;
    nd = 1e5;     % number of (accepted) grid search draws
    nBoot = 1000; % Number of bootstrap/projeciton draws
    nMC = 2000;   % Number of MC draws for coverage frequency
    nblocks = 50; % Number of workers on clusters
    level = 0.68; % confidence level of bounds on bounds 
    mineig = 0.001; % lower bound on smallest eigenvalue of sigma
    maxmod = 0.990; % upper bound on largest modulus of eigenvalue of companion matrix of A
 
    end
 
end

