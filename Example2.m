
UMPbaselineConfiguration = configFile;
UMPbaselineConfiguration.cum ='cum';
UMPbaselineConfiguration.SVARlabel = 'UMPbaseline';

load(['data' filesep 'UMP']) % the data is prepared outside this m-file

UMPlabelsOfTimeSeries = nams; % variable name tags
UMPtsInColumns = nums;

dataset = multivariateTimeSeries( UMPtsInColumns, UMPlabelsOfTimeSeries );
varEstimates = estimatedVecAR( UMPbaselineConfiguration, dataset);

restMat   = [ 1   0   1   0   1;
              2   0   1   0   1; 
              3   0  -1   0   1; 
              4   0   0   0   1]; 
IDscheme = IDassumptions( restMat); 

          
umpSVAR = SVAR( varEstimates, IDscheme);
%IRFCS = twoSidedIRFCS(umpSVAR,0.68,'analytic');
%plotPanel(IRFCS,'figures/UMP/CS')

%IRFCS.plotPanel()

% IRFid = umpSVAR.analytic.identifiedSet;
% IRFid.plotPanel

% 
% myexp = mcExperiment(umpSVAR);
% myexp.testCoverageTwoSidedIRFCSAnalytic(0.67)
% myexp.coverageTwoSidedIRFCSAnalytic.plotPanel('figures/UMP/');