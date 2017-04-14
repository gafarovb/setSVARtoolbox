classdef stochasticInequalities
    %stochasticInequalities Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        samples;
        originalSample;
        unitSphereGrid;
        config;
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
  
        
        function computeCSObjectiveFunctions(obj,gridPoint)
            %use bootstrap
        end
        function computeBonferroniCS(obj,level)
            alpha = 1-level;
            step1Level = 1 - obj.config.bonferroniStep1 * alpha;
            step2Level = 1 - (1 - obj.config.bonferroniStep1) * alpha;
            nGridpoints = obj.getNgridPoints;
            
            for i = 1: nGridpoints
                obj.testGridPoints( obj.unitSphereGrid(i,:)', step1Level)
                
%                 obj.computeCSObjectiveFunctions(gridPoint(i))
            end
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
        function resampledResidual = resampleResiduals(obj,gridPoint)
            Residual = obj.originalSample.computeEqualityResidualAtSphericalGridPoint(gridPoint);
            
            nSamples = size(obj.samples,2);
            nEqualities = size(Residual,1);

            resampledResidual = bootstrapSample(zeros(nEqualities,nSamples));
            independentNoise = obj.config.noiseStdToAvoidDeterministicConstraints * randn(nEqualities,nSamples);

            for i = 1: nSamples
                resampledResidual.values(:,i) =  obj.samples(i).computeEqualityResidualAtSphericalGridPoint(gridPoint) + independentNoise(:,i)  ;
            end
        end
        
        
        
        function resampledSlack = resampleSlacks(obj,gridPoint)
            Slack    = obj.originalSample.computeInequlaitySlackAtSphericalGridPoint(gridPoint);
            
            nSamples = size(obj.samples,2);
            nInequalities = size(Slack,1);
            
            resampledSlack =    bootstrapSample(zeros(nInequalities,nSamples));
            independentNoise = obj.config.noiseStdToAvoidDeterministicConstraints * randn(nInequalities,nSamples);

            for i = 1: nSamples
                resampledSlack.values(:,i)    =  obj.samples(i).computeInequlaitySlackAtSphericalGridPoint(gridPoint)  + independentNoise(:,i);
            end
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

