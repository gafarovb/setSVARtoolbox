classdef SVAR < VecAR
    
%%        SVAR class describes a set-identified SVAR model and offers tools
%       to construct point estimates and conduct inteference on the impulse 
%       response functions (IRF).
%
%       This version: 3.9.2017
%       Last edited by Bulat Gafarov 
%%    
    properties (Access = public) 
      label = 'Unknown' ; %% Model label, e.g. MSG 
      ID   =[] ; %% A class with restrictions
    end
    
    
%% *****************************************************************  
%  *****************************************************************  
%  *****************************************************************  
%  *****************************************************************      
%% *****************************************************************  

    methods
      function obj = SVAR(label,nLags)

         %  INPUTS: 
         %  label is the name of the model
         %  nLags is the number of lags in the model
      
         obj@VecAR(label,nLags); % create a reduced form VAR model
         obj.label = label ;
         obj.ID = restrictions(label); % read restrictions from a file 
      end
      
      function tempSuperClass(obj)
          
      end
      
      
      function n = nTS(obj)
         % Funciton nTS returns the number of time series
          n = obj.n;
      end 
      
    end
    
end

