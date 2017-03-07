classdef restrictions
    
    
%%      restrictions class describes restrictions on IRF in SVARS
%
%       This version: 3.7.2017
%       Last edited by Bulat Gafarov 
%%     


    properties
        restMatTMP;  % Matrix with columns: Var Hor S/Z Cum Shk
    end
    
    methods
        function obj = restrictions(label) 
           % read the restricitons matrix from a file  'restMat.dat'
            obj.restMatTMP = load([label filesep 'restMat.dat']);
        end
        function matrix = restMat(obj) 
            % read the restricitons matrix
            matrix = obj.restMatTMP;
        end
    end
    
end

