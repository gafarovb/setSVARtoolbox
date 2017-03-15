classdef IRF
    %IRF Summary of this class goes here
    %   Detailed explanation goes here
    
    
    properties (Access = public)
        Val = [];
        labelTS = 'Unnamed TS';
        label = 'Unknown IRF' ;
    end
    
    properties (Access = protected)
        MaxHorizons = 23;
    end
    
    methods
        function obj=IRF(labelTS,label,MaxHorizons)
            % inputs 
            %   labelTS      : name of ts corresponding to this IRF
            %   label       : name of the IRF. e.g. Upper bound delta CS
            %   MaxHorizons : max number of horizons 
            
            obj.labelTS = labelTS;
            obj.label   = label;
            obj.Val = zeros(1,MaxHorizons);
        end
    end
    
end

