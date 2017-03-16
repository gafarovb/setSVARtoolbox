classdef IRF
    %IRF Summary of this class goes here
    %   Detailed explanation goes here
    
    
    properties (Access = public)
        Val = [];
        labelTS = 'Unnamed TS';
        label = 'Unknown IRF' ;
    end
 
    
    methods
        function obj = IRF(vectorIRF, labelTS,label)
            % inputs
            %   vectorIRF   : IRF in form of a vector
            %   labelTS     : name of ts corresponding to this IRF
            %   label       : name of the IRF. e.g. Upper bound delta CS
          if nargin~=0  
            obj.labelTS = labelTS;
            obj.label   = label;
            obj.Val     = vectorIRF;
          end
        end
        function d = double(obj)
            d = obj.Val;
        end
        function obj = plus(obj1,obj2)
            a =  double(obj1);
            b =  double(obj2);
            obj = obj1;
            obj.Val = a + b;
        end
        function obj = minus(obj1,obj2)
            a =  double(obj1);
            b =  double(obj2);
            obj = obj1;
            obj.Val = a - b;
        end
        function nHorizons = MaxHorizons(obj)
            nHorizons = size(obj.Val,2)-1;
        end
        % plus
        % < >
    end
    
end

