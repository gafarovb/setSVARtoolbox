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

restMat   = [ 1   0   1   0   1;
              2   0   1   0   1; 
              3   0  -1   0   1; 
              4   0   0   0   1]; 
IDscheme = IDassumptions( restMat, shockName); 
          
level = 0.68;
          
umpSVAR = SVAR( varEstimates, IDscheme);

%%

load('./outputToolbox1/test_correct.mat')

idsetOldU = IRFcollection( Bounds.idU, UMPTSDiscription, []);
idsetOldL = IRFcollection( Bounds.idL, UMPTSDiscription, []);

oldIRF = join(idsetOldU, idsetOldL);

compareIRFIDSET = join(umpSVAR.analytic.identifiedSet , oldIRF);

idsetOldCSU = IRFcollection(Bounds.dltU,UMPTSDiscription,[]);
idsetOldCSL = IRFcollection(Bounds.dltL,UMPTSDiscription,[]);

oldIRF = join(idsetOldCSU, idsetOldCSL);
compareIRFCS = join(umpSVAR.analytic.twoSidedIRFCS(level), oldIRF);
 
plotPanel(join(compareIRFIDSET,compareIRFCS) ,['outputToolbox1' filesep])


 