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
      ID    ; %% A class with restrictions
      nLags ; %% The number of lags in the SVAR model
      names ; %% Labels for TS
    end
%%    
    properties (Access = protected)
      n = 0 ; % Number of time series
      T = 0 ;
      data ; 
    end
    
    
    methods
      function obj = SVAR(label,nLags)
         % This constructor opens a folder ./label/ and reads data
         % and identifying restrictions
         %  INPUTS: 
         %  label is the name of the model
         %  nLags is the number of lags in the model
        obj.label = label ;
        obj.ID = restrictions(label);
        if nargin<2
            obj.nLags = 12 ; 
        else
            obj.nLags = nLags ; 
        end
        % read data
        obj.data  = csvread([label filesep 'data.csv'],1);
        [obj.T,obj.n] = size( obj.data);   
        % read names of time series in a cell array
        fid       = fopen([label filesep 'data.csv']);
        obj.names = textscan(fid,[repmat('%[^,],',1,obj.n-1) '%[^,\r\n]'], 1);
        fclose(fid);
      end
      function n = nTS(obj)
         % Funciton nTS returns the number of time series
          n = obj.n;
      end 
      
    end
    
end

