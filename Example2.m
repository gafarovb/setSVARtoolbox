
UMPbaselineConfiguration = SVARconfiguration;
UMPbaselineConfiguration.isCumulativeIRF ='yes';
UMPbaselineConfiguration.SVARlabel = 'UMPbaseline';

load(['data' filesep 'UMP']) % the data is prepared outside this m-file

nTimeSeries = 4;
UMPlabelsOfTimeSeries = nams; % variable name tags
unitsOfMeasurement = repmat({'% change'},1,nTimeSeries);

UMPTSDiscription = [nams;unitsOfMeasurement];
UMPtsInColumns = nums;

dataset = multivariateTimeSeries( UMPtsInColumns, UMPTSDiscription );
varEstimates = estimatedVecAR( UMPbaselineConfiguration, dataset);

shockName = 'UMP';

%              TS  ho  sign cum
restMat   = [  1   0   1    0   ;
               2   0   1    0   ;
               3   0  -1    0   ;
               4   0  -1    0   ];
IDscheme = IDassumptions( restMat, shockName); 

level = 0.68;
          
umpSVAR = SVAR( varEstimates, IDscheme);

IRFidSet = estimatedIdentifiedSet(umpSVAR);

IRFCS = IRFtwoSidedCS(umpSVAR,level,'Analytic');


folderName = ['figures' filesep 'analytic' filesep ];
mkdir( folderName);
prefix = '';
pathAndPrefix_=[folderName prefix]; 
plotPanel([IRFidSet;IRFCS],pathAndPrefix_);



tic
IRFbonCS = IRFtwoSidedCS(umpSVAR, level,'Analytic');
toc


%IRFCS.plotPanel()

% IRFid = umpSVAR.analytic.identifiedSet;
% IRFid.plotPanel

% 
% myexp = mcExperiment(umpSVAR);
% myexp.testCoverageTwoSidedIRFCSAnalytic(0.67)
% myexp.coverageTwoSidedIRFCSAnalytic.plotPanel('figures/UMP/');