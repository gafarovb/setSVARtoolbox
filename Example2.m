%This version: May 3th, 2017.
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
               4   0   -1    0   ];
IDscheme = IDassumptions( restMat, shockName)
          
umpSVAR = SVAR( varEstimates, IDscheme)
 
IRFidSet = estimatedIdentifiedSet(umpSVAR);

IRFCS = IRFtwoSidedCS(umpSVAR,level,'Analytic');

folderName = ['figures' filesep 'analytic' filesep ];
mkdir( folderName);
prefix = '';
pathAndPrefix_=[folderName prefix]; 

panel1 = join(IRFidSet,IRFCS)
plotPanel( panel1,pathAndPrefix_);

tic
IRFbonCS = IRFtwoSidedCS(umpSVAR, level,'MI_Bonferroni');
toc

%% 
folderName = ['figures' filesep 'bonferroni' filesep ];
mkdir( folderName);

resultsFileName = [folderName 'autosave'];
save(resultsFileName,'IRFbonCS')
%%
prefix = '10e5_';
pathAndPrefix_=[folderName prefix]; 
IRFbonCS.setDescriptionField('type','--');
panel2 = join(IRFbonCS,IRFCS);

plotPanel( panel2,pathAndPrefix_);
% 
% myexp = mcExperiment(umpSVAR);
% myexp.testCoverageTwoSidedIRFCSAnalytic(0.67)
% myexp.coverageTwoSidedIRFCSAnalytic.plotPanel('figures/UMP/');