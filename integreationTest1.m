
UMPbaselineConfiguration = SVARconfiguration;
UMPbaselineConfiguration.isCumulativeIRF ='yes';
UMPbaselineConfiguration.SVARlabel = 'UMPbaseline';
load(['data' filesep 'UMP']) % the data is prepared outside this m-file


UMPlabelsOfTimeSeries = nams; % variable name tags
UMPtsInColumns = nums;

dataset = multivariateTimeSeries( UMPtsInColumns, UMPlabelsOfTimeSeries );
varEstimates = estimatedVecAR(UMPbaselineConfiguration, dataset);

restMat   = [ 1   0   1   0   1;
              2   0   1   0   1; 
              3   0  -1   0   1; 
              4   0   0   0   1]; 
IDscheme = IDassumptions( restMat); 
          
          
umpSVAR = SVAR( varEstimates, IDscheme);


%%

load('./outputToolbox1/test_correct.mat')

idsetOldU = IRFcollection(Bounds.idU,UMPlabelsOfTimeSeries,'Old point estimates upper bound');
idsetOldL = IRFcollection(Bounds.idL,UMPlabelsOfTimeSeries,'Old point estimates Lower bound');

oldIRF = [idsetOldU;idsetOldL];

oldIRF = oldIRF.setMarker('x'); 

compareIRF = [umpSVAR.analytic.identifiedSet;oldIRF];

compareIRF.plotPanel



idsetOldCSU = IRFcollection(Bounds.dltU,UMPlabelsOfTimeSeries,'Old point estimates upper bound');
idsetOldCSL = IRFcollection(Bounds.dltL,UMPlabelsOfTimeSeries,'Old point estimates Lower bound');

oldIRF = [idsetOldCSU;idsetOldCSL];

oldIRF = oldIRF.setMarker('x'); 

compareIRF = [umpSVAR.analytic.twoSidedIRFCS(0.68);oldIRF];

compareIRF.plotPanel('outputToolbox1')


 