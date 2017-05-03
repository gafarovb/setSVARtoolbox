classdef SVARconfiguration < handle
    %configFile contains specification for SVAR objects
    %   Detailed explanation goes here
    
    properties (Access = public )
        %% read specification
        assumptionsFilename = ['MSG' filesep 'restMat.dat'];
        dataFilenameCSV = ['MSG' filesep 'data.csv']; % T=165 for MSG
        SVARlabel = 'MSG';
        nLags = 12;
        nLagsMax = 24; % -    maximum possible number of lags in BIC/AIC/HQIC
        
        %% reduced var
        nNoncontemoraneousHorizons = 36;      % the number of horizons to compute the IRF
        scedasticity = 'homo'; % 'homo' = homoscedastic Omega, 'hetero' = heteroscedastic Omega
        babn1 = 1000;          % bab: number of bootstrap samples before bias correction
        babn2 = 1000;          % bab: number of bootstrap samples after  bias correction
   
        largeNumber = 100000;  % FIXME: constant c from the paper 
        smallNumber = 1e-6;
        %% Svar
        isCumulativeIRF = 'yes';  % 'yes' or 'no' 
        AndrewsSoaresKappa0 = 1.96;
        noiseStdToAvoidDeterministicConstraints = 1e-6;
        bonferroniStep1 = 0.5; 
        nGridPoints = 10000;
        nBootstrapSamples = 1000;
        
        %% Specify shock
        % solve for IRF bounds of specified shock
        % (by default, compute bounds for all shocks under gridsearch algorithm)
        masterSeed = 123456789; % Seed that generates seeds for every simulation
        MaxSimulations = 1000 ;  % number of Monte Carlo simulations
        shock = 1;

         
        
    end
    methods (Static)
        function tsReadyForAnalysis = prepareRawData(rawTSinColumns)
            tsReadyForAnalysis = MSG_demeanAndSumOfFirstAndFourthSeries(rawTSinColumns);
        end
        function f_T = andrewsSoaresTunningSequence(T)
            f_T = log(log(T))/sqrt(T);
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