classdef estimatedVecAR < VecAR
    %VecAR is a class for a time series collection
    %   Detailed explanation goes here
    
    %%
    properties (Access = protected)
        dataSample;
        estimates;
        
    end
    

    methods  % constructors
        function obj = estimatedVecAR( config,dataset)
            if nargin<1
                warning('Configruation  is not provided. Using the default VAR model and dataset')
                obj.config = configFile;
            else
                obj.config = config;
            end
            obj.nLags = obj.config.nLags ;
            obj.scedasticity = obj.config.scedasticity ;
            
            if nargin<2
                warning('Dataset is not provided. Using the default dataset')
                [tsInColumns, labelsOfTimeSeries ] = obj.readDataFromFile(obj.config);
                dataset = multivariateTimeSeries( tsInColumns, labelsOfTimeSeries );
            end
            
            obj.dataSample = dataset;
            obj.estimates = LSestimatesVAR( obj.dataSample, obj.nLags,obj.config.nNoncontemoraneousHorizons, obj.scedasticity);
        end
        
        function [tsInColumns, labelsOfTimeSeries ]= readDataFromFile(obj,config)
            tsInColumns  = csvread(config.dataFilenameCSV,1);
            tsInColumns  = config.prepareRawData( tsInColumns);
            nTS = size(tsInColumns,2);
            labelsOfTimeSeries = obj.readCSVheader( config.dataFilenameCSV, nTS);
        end
    end
    methods  % basic characteristics
        function namesOfTS = getNames(obj)
            namesOfTS = obj.dataSample.getNames;
        end
        function n = getN(obj)
            n = obj.dataSample.countTS;
        end
        function T = getT(obj)
            T = obj.dataSample.countTimePeriods;
        end
        function d = countParameters(obj)
            n = obj.getN ;
            p = obj.nLags;
            d = n * n * p + n * (n+1) / 2;
        end
    end
    methods  % estimates
        function Sigma = getSigma(obj)
            Sigma = obj.estimates.getSigma;
        end
        function thetaHat = getTheta(obj)
            thetaHat = obj.estimates.getThetaHat;
        end
        function OmegaT = getCovarianceOfThetaT(obj)
            OmegaT = obj.estimates.getOmega / obj.getT ;
        end
    end
    methods  % Represetnations
        function VMA_ts_sh_ho = getVMA_ts_sh_ho(obj)
            VMA_ts_sh_ho = obj.estimates.getVMA_ts_sh_ho( obj.config.nNoncontemoraneousHorizons);
        end
        function G = getVMADerivatives_ts_sh_ho_dAL(obj)
            G = obj.estimates.getVMADerivatives;
        end
    end
    methods % analysis of the VARs
        function stationarityTest(obj)
            % Checks whether reduced-form VAR model is covariance stationary
            VecAR.stationarity(obj.estimates.getAL)
        end
        function [swartzLags, akaikeLags,hannanLags] = optimalNlagsByIC(obj)
            % This function selects number of lags in the VAR model for "series"
            % based on AIC,BIC and HQ information criteria
            %
            % This version: March 23, 2017
            nLagsMax = obj.config.nLagsMax;
            aic = zeros(nLagsMax,1);
            bic = zeros(nLagsMax,1);
            hqic = zeros(nLagsMax,1);
            
            for p = 1:nLagsMax
                estimatesForPlags = LSestimatesVAR( obj.dataSample, p,obj.config.nNoncontemoraneousHorizons, obj.scedasticity);
                [bic(p), aic(p), hqic(p)] = getBicAicHQic(estimatesForPlags);
                obj.estimates = LSestimatesVAR( obj.dataSample, obj.nLags ,obj.config.nNoncontemoraneousHorizons, obj.scedasticity );
            end
            
            [~, akaikeLags] = min(aic);
            [~, swartzLags] = min(bic);
            [~, hannanLags] = min(hqic);
            
        end
    end
        methods (Static)
        function firstRow = readCSVheader(dataFilenameCSV,nColumns)
            fid       = fopen(dataFilenameCSV);
            firstRow = textscan(fid,[repmat('%[^,],',1,nColumns-1) '%[^,\r\n]'], 1);
            fclose(fid);
        end
    end
end




