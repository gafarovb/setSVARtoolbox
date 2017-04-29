
UMPbaselineConfiguration = configFile;
UMPbaselineConfiguration.cum ='cum';
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
         
          


