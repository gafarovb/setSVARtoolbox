classdef restrictions
    
    
%%      restrictions class describes restrictions on IRF in SVARS
%
%       This version: 3.7.2017
%       Last edited by Bulat Gafarov 
%%     


    properties
        restMat;  % Matrix with columns: Var Hor S/Z Cum Shk
    end
    
    methods
        function obj = restrictions(label) 
            obj.restMat = load([label filesep 'restrictions.dat']);
        end
    end
    
end

