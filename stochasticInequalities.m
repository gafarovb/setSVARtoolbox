classdef stochasticInequalities < handle
    %stochasticInequalities Summary of this class goes here
    %   Detailed explanation goes here
    properties (Access = private)
        samples;
        originalSample;
        unitSphereGrid;
        config;
        waitBarWindow;
    end
    
    methods
        function obj = stochasticInequalities(SVARobj)
            obj.config = SVARobj.getConfig;
            obj.originalSample = SVARobj;
            obj.samples = SVARobj.generateSamplesFromAsymptoticDistribution( obj.getNsamples);
            obj.unitSphereGrid = obj.generateNdimUnitSphereGrid(SVARobj.getN, obj.getNgridPoints);
            
            obj.waitBarWindow = waitBarCustomized(  obj.getNgridPoints);
            obj.waitBarWindow.setMessage('Computing Bonferroni CS, stage 1 ');

        end
        function b = getNsamples(obj)
            b = obj.config.nBootstrapSamples;
        end
        function g = getNgridPoints(obj)
            g = obj.config.nGridPoints;
        end
    end
    methods
        function [lowerBoundIRF, upperBoundsIRF] = computeBonferroniCS(obj,level)
            alpha = 1-level;
            step1Level = 1 - obj.config.bonferroniStep1 * alpha;
            step2Level = 1 - (1 - obj.config.bonferroniStep1) * alpha;
            nGridpoints = obj.getNgridPoints;
            
            gridPointInCS = true(nGridpoints,1);
            
            for i = 1: nGridpoints
                gridPointInCS(i) = obj.testGridPoints( obj.unitSphereGrid(i,:)', step1Level);
                obj.waitBarWindow.showProgress(i);
            end
            
            for i = nGridpoints: -1: 1
                if  gridPointInCS(i)
                    [lowerBounds(:,:,i), upperBounds(:,:,i) ] = obj.computeCSObjectiveFunctions( obj.unitSphereGrid(i,:)', step2Level);
                end
            end
            
            namesTS = obj.originalSample.getNamesOfTS;
            label = ['Bonferroni 2nd step CS with p= ' num2str(level) ];
            lowerBoundIRF = IRFcollection( nanmin(lowerBounds,[],3), namesTS, label );
            upperBoundsIRF = IRFcollection( nanmax(upperBounds,[],3), namesTS, label );
        end
        function isInTheCS = testGridPoints(obj,gridPoint,level)
            resampledResidual = obj.resampleResiduals( gridPoint);
            resampledSlack = obj.resampleSlacks( gridPoint);
            
            Slack    = obj.originalSample.computeInequlaitySlackAtSphericalGridPoint(gridPoint);
            SlackStandardized = Slack./std(resampledSlack);
            activeSet = obj.generalizedMomentSelection(SlackStandardized);
            
            Residual = obj.originalSample.computeEqualityResidualAtSphericalGridPoint(gridPoint);
            ResidualStandardized = Residual./std(resampledResidual);
            mmm = obj.modifiedMethodOfMoments( ResidualStandardized, SlackStandardized(activeSet));
            
            resampledSlackActive = bootstrapSample(resampledSlack.values(activeSet,:));
            quantileMMM = obj.computeQuantileOfModifiedMethodOfMoments(resampledResidual,resampledSlackActive, level)  ;
            
            isInTheCS = (mmm <= quantileMMM);
        end
        function activeSet = generalizedMomentSelection(obj,standardizefSlack)
            T = obj.originalSample.getT;
            config = obj.originalSample.getConfig;
            f_T = config.AndrewsSoaresKappa0 * config.andrewsSoaresTunningSequence(T);
            slackSet = (standardizefSlack>=f_T);
            activeSet = ~slackSet;
        end
        function quantileMMM = computeQuantileOfModifiedMethodOfMoments(obj,resids,slacks,level)
            mmmResampled = bootstrapSample(obj.modifiedMethodOfMoments(resids.studentized,slacks.studentized));
            quantileMMM = mmmResampled.quantile(level);
        end
        function [lowerBound, upperBound] = computeCSObjectiveFunctions(obj,gridPoint, level)
            lowerQ = (1 - level)/2;
            upperQ = 1 - lowerQ;
            resampledObjectiveFunctions = obj.resampleObjectiveFunctions( gridPoint);

            lowerBound = resampledObjectiveFunctions.quantile(lowerQ);
            upperBound = resampledObjectiveFunctions.quantile(upperQ);
        end
    end
    methods
        function resampledObjectiveFunctions = resampleObjectiveFunctions(obj,gridPoint)
            nSamples = obj.getNsamples;
            
            originalObjectiveFunctions = obj.originalSample.computeObjectiveFunctionsAtSphericalGridPoint( gridPoint);
            shape = size(originalObjectiveFunctions);
            
            samplingDimension = 3;
            resampledObjectiveFunctions = bootstrapSample( zeros(shape(1),shape(2), nSamples), samplingDimension);
            
            for i = 1: nSamples
                resampledObjectiveFunctions.values(:,:,i) =  obj.samples(i).computeObjectiveFunctionsAtSphericalGridPoint(gridPoint)  ;
            end
            
        end
        function resampledResidual = resampleResiduals(obj,gridPoint)
            nSamples = obj.getNsamples;
            
            Residual = obj.originalSample.computeEqualityResidualAtSphericalGridPoint(gridPoint);
            nEqualities = size(Residual,1);
            
            resampledResidual = bootstrapSample(zeros(nEqualities,nSamples));
            independentNoise = obj.config.noiseStdToAvoidDeterministicConstraints * randn(nEqualities,nSamples);
            
            for i = 1: nSamples
                resampledResidual.values(:,i) =  obj.samples(i).computeEqualityResidualAtSphericalGridPoint(gridPoint) + independentNoise(:,i)  ;
            end
        end
        function resampledSlack = resampleSlacks(obj,gridPoint)
            
            nSamples =  obj.getNsamples;
            
            Slack    = obj.originalSample.computeInequlaitySlackAtSphericalGridPoint(gridPoint);
            nInequalities = size(Slack,1);
            
            resampledSlack =    bootstrapSample( zeros(nInequalities, nSamples));
            independentNoise = obj.config.noiseStdToAvoidDeterministicConstraints * randn(nInequalities,nSamples);
            
            for i = 1: nSamples
                resampledSlack.values(:,i)    =  obj.samples(i).computeInequlaitySlackAtSphericalGridPoint(gridPoint)  + independentNoise(:,i);
            end
        end
    end
    methods(Static)
        function gridOnSphere = generateNdimUnitSphereGrid(dimension,nGridPoints )
            gridOnSphere = randn(nGridPoints,dimension);
            gridNormColumn = sqrt(   sum(  gridOnSphere.*gridOnSphere,2));
            gridNorm = repmat(gridNormColumn,1,dimension);
            gridOnSphere = gridOnSphere./gridNorm;
        end
        function mmm = modifiedMethodOfMoments(resid,slack)
            jointVector = [resid;slack .* (slack<=0)];
            mmm = sum(jointVector.*jointVector,1);
        end
    end
end

