classdef classMSG < SVAR  
    %CLASSMSG Summary of this class goes here
    %   Detailed explanation goes here
    
    
    methods 
        function obj = classMSG(nLags)
        obj@SVAR('MSG',nLags);
           
        %% adjustments
        
        obj.data(:,4) = obj.data(:,4)+obj.data(:,1);
        obj.data = obj.data - ones(obj.T,1)*mean(obj.data);  % demean series 

        end
    end
    
end

