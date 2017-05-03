classdef stochasticInequalities 
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
        end
        function b = getNsamples(obj)
            b = obj.config.nBootstrapSamples;
        end
        function g = getNgridPoints(obj)
            g = obj.config.nGridPoints;
        end
    end
    methods
        function [lowerBoundIRF,upperBoundsIRF ] = computeBonferroniCS(obj,level)
            alpha = 1-level;
            
            step1Level = 1 - obj.config.bonferroniStep1 * alpha;
            step2Level = 1 - (1 - obj.config.bonferroniStep1) * alpha;
            
            gridPointInCS = obj.testAllPoints(step1Level);
            obj.assertIfNoFeasiblePointsFound(gridPointInCS);
            
            [lowerBounds, upperBounds ] = obj.computeConditionalIRFCS( gridPointInCS,step2Level);
            
            IRFDescription = descriptionBonferroniCS(obj,level);
            TSDescription = obj.originalSample.getTSDescription;
            lowerBoundIRF = IRFcollection( nanmin(lowerBounds(:,:,gridPointInCS),[],3), TSDescription, IRFDescription );
            upperBoundsIRF = IRFcollection( nanmax(upperBounds(:,:,gridPointInCS),[],3),TSDescription, IRFDescription );
        end
        
        function IRFDescription = descriptionBonferroniCS(obj,level)
            levelStr =   num2str(level);
            SVARobj = obj.originalSample;
            
            IRFDescription = struct('tag','noTag',... % can be used in filenames
                'legend',['Bonferroni 2nd step CS with p= ' levelStr ],...
                'shock', SVARobj.getShockLabel,...
                'SVARmodel', SVARobj.label,...
                'type' , 'CS');
            
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
        function isInTheCS = testGridPoint(obj,gridPoint,level)
            resampledResidual = obj.resampleResiduals( gridPoint);
            resampledSlack = obj.resampleSlacks( gridPoint);
            
            Slack    = obj.originalSample.MIframework.computeInequlaitySlackAtSphericalGridPoint(gridPoint);
            SlackStandardized = Slack./std(resampledSlack);
            activeSet = obj.generalizedMomentSelection(SlackStandardized);
            
            Residual = obj.originalSample.MIframework.computeEqualityResidualAtSphericalGridPoint(gridPoint);
            ResidualStandardized = Residual./std(resampledResidual);
            mmm = obj.modifiedMethodOfMoments( ResidualStandardized, SlackStandardized(activeSet));
            
            resampledSlackActive = bootstrapSample(resampledSlack.values(activeSet,:));
            quantileMMM = obj.computeQuantileOfModifiedMethodOfMoments(resampledResidual,resampledSlackActive, level)  ;
            
            isInTheCS = (mmm <= quantileMMM);
        end
        function gridPointInCS = testAllPoints(obj,step1Level)
            nGridpoints = obj.getNgridPoints;
            obj.waitBarWindow = waitBarCustomized(  obj.getNgridPoints);
            gridPointInCS = false(nGridpoints,1);
            
            if isempty(gcp('nocreate'))
                obj.waitBarWindow.setMessage('Computing Bonferroni CS, Step 1');
                for i = 1:nGridpoints
                    gridPointInCS(i) = obj.testGridPoint( obj.unitSphereGrid(i,:)', step1Level);
                    obj.waitBarWindow.showProgress( i);
                end
            else
                obj.waitBarWindow.setMessage('Computing Bonferroni CS, Step 1, parallel mode');
                obj.waitBarWindow.showProgress( 1);
                parfor i = 1:nGridpoints
                    gridPointInCS(i) = testGridPoint(obj, obj.unitSphereGrid(i,:)', step1Level);
                end
            end
            delete(obj.waitBarWindow);
        end
        function [lowerBounds, upperBounds ] = computeConditionalIRFCS(obj,gridPointInCS,step2Level)
            lowerBounds = obj.emptyIRF;
            upperBounds = obj.emptyIRF;
            nGridpoints = size(gridPointInCS,1);
            
            iFeasiblePoint = 0;
            obj.waitBarWindow = waitBarCustomized( sum(gridPointInCS));
            if isempty(gcp('nocreate'))
                obj.waitBarWindow.setMessage(['Computing Bonferroni CS, Step 2 (' num2str(sum(gridPointInCS)) ' feasible point found)']);
                for i = 1: nGridpoints
                    if gridPointInCS(i)
                        [lowerBounds(:,:,i), upperBounds(:,:,i) ] = obj.computeCSObjectiveFunctionsAtGridPoint( obj.unitSphereGrid(i,:)', step2Level);
                        iFeasiblePoint = iFeasiblePoint+1;
                        obj.waitBarWindow.showProgress( iFeasiblePoint);
                    end
                end
            else
                obj.waitBarWindow.setMessage(['Computing Bonferroni CS, Step 2 (' num2str(sum(gridPointInCS)) ' feasible point found), parallel mode']);
                obj.waitBarWindow.showProgress( 1);

                parfor i = 1:nGridpoints
                    [lowerBounds(:,:,i), upperBounds(:,:,i) ] = computeCSObjectiveFunctionsAtGridPoint(obj, obj.unitSphereGrid(i,:)', step2Level);
                end
            end
            delete(obj.waitBarWindow);
        end
        function [lowerBound, upperBound] = computeCSObjectiveFunctionsAtGridPoint(obj,gridPoint, level)
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
            
            originalObjectiveFunctions = obj.originalSample.MIframework.computeObjectiveFunctionsAtSphericalGridPoint( gridPoint);
            shape = size(originalObjectiveFunctions);
            
            samplingDimension = 3;
            resampledObjectiveFunctions = bootstrapSample( zeros(shape(1),shape(2), nSamples), samplingDimension);
            
            for i = 1: nSamples
                resampledObjectiveFunctions.values(:,:,i) =  obj.samples(i).MIframework.computeObjectiveFunctionsAtSphericalGridPoint(gridPoint)  ;
            end
            
        end
        function resampledResidual = resampleResiduals(obj,gridPoint)
            nSamples = obj.getNsamples;
            
            Residual = obj.originalSample.MIframework.computeEqualityResidualAtSphericalGridPoint(gridPoint);
            nEqualities = size(Residual,1);
            
            resampledResidual = bootstrapSample(zeros(nEqualities,nSamples));
            independentNoise = obj.config.noiseStdToAvoidDeterministicConstraints * randn(nEqualities,nSamples);
            
            for i = 1: nSamples
                resampledResidual.values(:,i) =  obj.samples(i).MIframework.computeEqualityResidualAtSphericalGridPoint(gridPoint) + independentNoise(:,i)  ;
            end
        end
        function resampledSlack = resampleSlacks(obj,gridPoint)
            
            nSamples =  obj.getNsamples;
            
            Slack    = obj.originalSample.MIframework.computeInequlaitySlackAtSphericalGridPoint(gridPoint);
            nInequalities = size(Slack,1);
            
            resampledSlack =    bootstrapSample( zeros(nInequalities, nSamples));
            independentNoise = obj.config.noiseStdToAvoidDeterministicConstraints * randn(nInequalities,nSamples);
            
            for i = 1: nSamples
                resampledSlack.values(:,i)    =  obj.samples(i).MIframework.computeInequlaitySlackAtSphericalGridPoint(gridPoint)  + independentNoise(:,i);
            end
        end
        function emptyArray =  emptyIRF(obj)
            emptyArray =  NaN(obj.originalSample.getN, obj.originalSample.getMaxHorizons, obj.getNgridPoints);
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
        function assertIfNoFeasiblePointsFound(gridPointInCS)
            if sum(gridPointInCS)==0
                error('SVAR_CSBonferroni:noFeasiblePointFound','No feasible point found. Try a larger number of the gridpoints.');
            end
        end
       
    end
end

