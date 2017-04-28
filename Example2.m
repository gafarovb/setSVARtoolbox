
UMPbaselineConfiguration = configFile;
UMPbaselineConfiguration.cum ='cum';
UMPbaseline.SVARlabel = 'UMPbaseline';
load(['data' filesep 'UMP']) % the data is prepared outside this m-file


UMPlabelsOfTimeSeries = nams; % variable name tags
UMPtsInColumns = nums;

varEstimates = estimatedVecAR(UMPbaselineConfiguration, UMPtsInColumns, UMPlabelsOfTimeSeries);

restMat   = [ 1   0   1   1   1;
              2   0   1   1   1; 
              3   0  -1   1   1; 
              4   0   0   1   1]; 
          
umpSVAR = SVAR(varEstimates, restMat);
         
          


