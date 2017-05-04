
UMPbaselineConfiguration = SVARconfiguration;
UMPbaselineConfiguration.nGridPoints = 100000;
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

 
tic
IRFbonCS = IRFtwoSidedCS(umpSVAR, level,'MI_Bonferroni');
toc

%% 
folderName = ['figures' filesep 'bonferroni' filesep ];
mkdir( folderName);
prefix = ['10e' num2str( log10(UMPbaselineConfiguration.nGridPoints)) '_'];
resultsFileName = [folderName prefix 'autosave'];
save(resultsFileName,'IRFbonCS')
%%

pathAndPrefix_=[folderName prefix]; 
IRFbonCS = IRFbonCS.setDescriptionField('type','+');
panel2 = join(IRFbonCS,IRFCS);

plotPanel( panel2,pathAndPrefix_);
% 
% myexp = mcExperiment(umpSVAR);
% myexp.testCoverageTwoSidedIRFCSAnalytic(0.67)
% myexp.coverageTwoSidedIRFCSAnalytic.plotPanel('figures/UMP/');