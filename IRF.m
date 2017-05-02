classdef IRF 
    %IRF Summary of this class goes here
    %   Detailed explanation goes here
    
    
    properties (Access = public)
        Values = [];
        labelTS = 'Unnamed TS';
        label = 'Unknown IRF' ;
        model = 'UnknownModel';
        unitsOfMeasurement = '% change';
        markerString = '-';
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
            obj.Values   = vectorIRF;
          end
        end
        function d = double(obj)
            d = obj.Values;
        end
        function obj = plus(obj1,obj2)
            [obj,a,b] = genericOperation(obj1,obj2);
             obj = obj.setValues(a + b);
        end
        function obj = mtimes(obj1,obj2)
            [obj,a,b] = genericOperation(obj1,obj2);
            obj = obj.setValues(a.*b);
        end
        function obj = minus(obj1,obj2)
            [obj,a,b] = genericOperation(obj1,obj2);
            obj = obj.setValues(a - b);
        end
        function obj = le(obj1,obj2)
            [obj,a,b] = genericOperation(obj1,obj2);
            obj = obj.setValues(logical(a <= b));
        end
        function obj = ge(obj1,obj2)
            [obj,a,b] = genericOperation(obj1,obj2);
            obj = obj.setValues(logical(a >= b));
        end
        function obj = setValues(obj,values)
            obj.Values = values;
        end
        function nHorizons = nNoncontemoraneousHorizons(obj)
            nHorizons = size(obj.Values,2)-1;
        end
        function obj = setLabel(obj,label)
            obj.label = label;
        end 
        function obj = setUnits(obj,unitString)
            obj.unitsOfMeasurement = unitString;
        end 
        function fileName = getLabelTS(obj)
            fileName =  char(obj.labelTS{1})  ;
        end
  
        function figHandle = plot(obj)
            
            maxHorizon = obj.nNoncontemoraneousHorizons;
            figHandle = plot(0:1:maxHorizon, obj.Values,['k'  obj.markerString ],'LineWidth',2);
            title(obj.labelTS{1},'Interpreter','tex','FontSize',12); 
            hold on; 
            ylabel(obj.unitsOfMeasurement); 
            xlabel('Months after shock');
            set(gca,'LineWidth',2.0);
            grid on; 
            box off; 
        end
    end
    
end

function [obj,a,b] = genericOperation(obj1,obj2)
            a =  double(obj1);
            b =  double(obj2);
            obj = obj2;
end
