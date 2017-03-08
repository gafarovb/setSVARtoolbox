classdef VecAR 
    %VAR is a class for a time series collection 
    %   Detailed explanation goes here
 
    properties (Access = public) 
      nLags = 0 ; %% The number of lags in the SVAR model
      names = []; %% Labels for TS
    end
  
%%  

    properties (Access = protected)
      n = 0 ; % Number of time series
      T = 0 ;
      data = [] ; 
    end
    
    
    methods
        function obj = VecAR(label,nLags)
         % This constructor opens a folder ./label/ and reads data
         % and identifying restrictions    
            
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
        
    end
    
end

