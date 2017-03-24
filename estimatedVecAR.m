classdef estimatedVecAR < VecAR
    %VecAR is a class for a time series collection
    %   Detailed explanation goes here

    %%
    properties (Access = protected)
        dataSample; 
        estimates;
    end
    
    methods
        function obj = estimatedVecAR
            obj.config = configFile;
            obj.nLags = obj.config.nLags ;
            obj.dataSample = multivariateTimeSeries;
            obj.estimates = LSestimatesVAR( obj.dataSample, obj.nLags);
        end
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
                estimatesForPlags = LSestimatesVAR( obj.dataSample, p);
                [bic(p), aic(p), hqic(p)] = getBicAicHQic(estimatesForPlags);
                obj.estimates = LSestimatesVAR( obj.dataSample, obj.nLags);
            end
            
            [~, akaikeLags] = min(aic);
            [~, swartzLags] = min(bic);
            [~, hannanLags] = min(hqic);
            
        end
        function n = getN(obj)
            n = obj.dataSample.countTS;
        end
        function d = countParameters(obj)
            n = obj.getN ;
            p = obj.nLags;
            d = n * n * p + n * (n+1) / 2;
        end
        function thetaCov = getOmega(obj)
            thetaCov = obj.estimates.getOmega;
        end
        function thetaHat = getTheta(obj)
            thetaHat = obj.estimates.getThetaHat;
        end
        function VMA_ts_sh_ho = getVMA_ts_sh_ho(obj)
            VMA_ts_sh_ho = obj.estimates.getVMA_ts_sh_ho;
        end
        function G = getVMADerivatives(obj)
            G = obj.estimates.getVMADerivatives;
        end
        function T = getT(obj)
            T = obj.dataSample.countTimePeriods;
        end
    end
    
end




