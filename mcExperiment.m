classdef mcExperiment < handle
    %MCSVAR Summary of this class goes here
    % This class describes a Monte Carlo exmperiment
    %   Detailed explanation goes here
    
    properties
        design;                  % SVAR object with design
        coverageOnesidedUpperIRFCSAnalytic;
        coverageOnesidedLowerIRFCSAnalytic;
        coverageTwoSidedIRFCSAnalytic;
    end
    
    properties (Access = private)
        config;
        samples;                 % Array of SVAR objects with generated samples
        waitBarWindow;
    end
    
    methods
        function obj = mcExperiment(SVARobj)
            obj.design = SVARobj;
            obj.config = SVARobj.getConfig;
            obj.waitBarWindow = waitBarCustomized(  obj.getNumberOfSimulations);
            obj.samples = SVARobj.generateSamplesFromAsymptoticDistribution( obj.getNumberOfSimulations);
        end     
        function MaxSimulations = getNumberOfSimulations(obj)
            MaxSimulations = obj.config.MaxSimulations;
        end

        function  testCoverageOnesidedUpperIRFCSAnalytic(obj,nominalLevel)
            H0_upperBounds = obj.design.analytic.onesidedUpperIRFHat;
            nShocks = size(H0_upperBounds,2);
            H0isNotRejected(obj.getNumberOfSimulations,nShocks) = IRFcollection;
            for i = 1: obj.getNumberOfSimulations
                H0isNotRejected(i,:) =  (obj.samples(i).analytic.onesidedUpperIRFCS(nominalLevel) >= H0_upperBounds);
                obj.MCwaitBar(i);
            end
            obj.coverageOnesidedUpperIRFCSAnalytic = obj.computeCoverageFrequency(H0isNotRejected);
            obj.coverageOnesidedUpperIRFCSAnalytic = obj.coverageOnesidedUpperIRFCSAnalytic.setMarker('o');
            obj.coverageOnesidedUpperIRFCSAnalytic = obj.coverageOnesidedUpperIRFCSAnalytic.setLabel(['analytic upper one-sided CS cowerage with nominal p=',num2str(nominalLevel)]);
        end
        function  testCoverageOnesidedLowerIRFCSAnalytic(obj,nominalLevel)
            H0_lowerBounds = obj.design.analytic.onesidedLowerIRFHat;
            nShocks = size(H0_lowerBounds,2);
            H0isNotRejected(obj.getNumberOfSimulations,nShocks) = IRFcollection;
            for i = 1: obj.getNumberOfSimulations
                H0isNotRejected(i,:) =  (obj.samples(i).analytic.onesidedLowerIRFCS(nominalLevel) <= H0_lowerBounds);
                obj.MCwaitBar(i);
            end
            obj.coverageOnesidedLowerIRFCSAnalytic =  obj.computeCoverageFrequency(H0isNotRejected);
            obj.coverageOnesidedLowerIRFCSAnalytic = obj.coverageOnesidedLowerIRFCSAnalytic.setMarker('o');
            obj.coverageOnesidedLowerIRFCSAnalytic = obj.coverageOnesidedLowerIRFCSAnalytic.setLabel(['analytic lower one-sided CS cowerage with nominal p=',num2str(nominalLevel)]);
        end      
        function  testCoverageTwoSidedIRFCSAnalytic(obj,nominalLevel)
            H0_upperBounds = obj.design.analytic.onesidedUpperIRFHat;
            H0_lowerBounds = obj.design.analytic.onesidedLowerIRFHat;
            nShocks = size(H0_upperBounds,2);
            H0isNotRejected(obj.getNumberOfSimulations,nShocks) = IRFcollection;
            twoSidedLevel = 1 - (1 - nominalLevel)/2;
            
            for i = 1: obj.getNumberOfSimulations
                H0isNotRejected(i,:) =  (obj.samples(i).analytic.onesidedUpperIRFCS(twoSidedLevel) >= H0_upperBounds)* ...
                                        (obj.samples(i).analytic.onesidedLowerIRFCS(twoSidedLevel) <= H0_lowerBounds);
                obj.MCwaitBar(i);
            end
            obj.coverageTwoSidedIRFCSAnalytic =  obj.computeCoverageFrequency(H0isNotRejected);
            obj.coverageTwoSidedIRFCSAnalytic = obj.coverageTwoSidedIRFCSAnalytic.setMarker('o');
            obj.coverageTwoSidedIRFCSAnalytic = obj.coverageTwoSidedIRFCSAnalytic.setLabel( ['analytic two-sided CS cowerage with nominal p=',num2str(nominalLevel)]);
        end
        function  coverageFreqIRF = computeCoverageFrequency(obj,isCoveredIRFCollection)
            coverageFreqIRF =  isCoveredIRFCollection(1,:);
            for i = 2: obj.getNumberOfSimulations
                coverageFreqIRF = isCoveredIRFCollection(i,:) + coverageFreqIRF;
            end
            coverageFreqIRF =  (1 / obj.getNumberOfSimulations) * coverageFreqIRF; % don't change the order of multiplicaiton
        end
        function  MCwaitBar(obj,step)
            obj.waitBarWindow.showProgress(step);
        end
    end
    
end

