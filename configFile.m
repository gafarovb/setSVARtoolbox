classdef configFile < handle
    %configFile contains specification for SVAR objects
    %   Detailed explanation goes here
    
    properties (Access = public,Constant )
        %% read specification
        assumptionsFilename = ['MSG' filesep 'restMat.dat'];
        dataFilenameCSV = ['MSG' filesep 'data.csv']; % T=165
        label = 'MSG';
        nLags = 12;
        nLagsMax = 24; % -    maximum possible number of lags in BIC/AIC/HQIC
        
        %% reduced var
        MaxHorizons = 23;      % the number of horizons to compute the IRF
        scedasticity = 'homo';      % 'homo' = homoscedastic Omega, 'hetero' = heteroscedastic Omega
        babn1 = 1000;          % bab: number of bootstrap samples before bias correction
        babn2 = 1000;          % bab: number of bootstrap samples after  bias correction
   
        largeNumber = 100000;  % FIXME: constant c from the paper 
        smallNumber = 1e-6;
        %% Svar
        cum = '';  % 'cum' = cumulative IRF responses

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
    methods (Static)
        function tsReadyForAnalysis = prepareRawData(rawTSinColumns)
            tsReadyForAnalysis = MSG_demeanAndSumOfFirstAndFourthSeries(rawTSinColumns);
        end
    end
        

end


function tsReadyForAnalysis  = MSG_demeanAndSumOfFirstAndFourthSeries(rawTSinColumns)

    tsReadyForAnalysis = rawTSinColumns;
    tsReadyForAnalysis(:,4) = tsReadyForAnalysis(:,4) + tsReadyForAnalysis(:,1);
    tsReadyForAnalysis = demeanColumns(tsReadyForAnalysis);  
    
end
function columnsDemeaned = demeanColumns(columnsWithMeans)
    T = size(columnsWithMeans,1);
    columnsDemeaned = columnsWithMeans - ones(T,1)*mean(columnsWithMeans);  
end