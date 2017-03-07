classdef SVAR
    
%%        SVAR class describes a set-identified SVAR model and offers tools
%       to construct point estimates and conduct inteference on the impulse 
%       response functions (IRF).
%
%       This version: 3.7.2017
%       Last edited by Bulat Gafarov 
%%    
    properties (Access = public) 
      label = 'Unknown' ; %% Model label, e.g. MSG 
      idRestrictions ; %% A class with restrictions
      nLags ; %% the number of lags in the SVAR model
    end
%%    
    properties (Access = protected)
      n = 1 ; % Number of time series     
    end
    
    
    methods
      function obj = SVAR(label,nLags)
         % This constructor opens a folder ./label/ and reads data and the 
         
         % label is the name of the model.
        obj.label = label ;
        obj.idRestrictions;
        obj.nLags = nLags ; 
      end
      function n = nTS(obj)
         % Funciton nTS returns the number of time series
          n = obj.n;
      end 
      
    end
    
end

