classdef IRF 
    %IRF Summary of this class goes here
    %   Detailed explanation goes here
    
    
    properties (Access = public)
        Values = [];
        TSdescription = ['Unnamed TS' 'Units'] ;
        description = struct('tag','noTag',... % can be used in filenames
                'legend','Unknown IRF',...
                'shock', 'Unknown shock',...
                'SVARmodel', 'Unkown model',...
                'type' , 'Unknown');
    end
 
    
    methods
        function obj = IRF(vectorIRF, TSdescription, description)
            % inputs
            %   vectorIRF   : IRF in form of a vector
            %   labelTS     : name of ts corresponding to this IRF
            %   description       
            
          if nargin~=0  
              if ~isempty(TSdescription)
                  obj.TSdescription = TSdescription;
              end
              if ~isempty(description)
                  obj.description   = description;
              end
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
        function obj = setDescription(obj,description)
            obj.label = description;
        end 
        function obj = setDescriptionField(obj,fieldName, newString)
            obj.description.(fieldName) = newString;
        end 
        function fileName = getLabelTS(obj)
            fileName =  char(obj.TSdescription{1,1})  ;
        end
        function figHandle = plot(obj)
            maxHorizon = obj.nNoncontemoraneousHorizons;
            figHandle = plot(0:1:maxHorizon, obj.Values,['k'  obj.getMarker ],'LineWidth',2);
            title(obj.TSdescription{1,1},'Interpreter','tex','FontSize',12); 
            hold on; 
            ylabel(obj.TSdescription{2,1}); 
            xlabel('Months after shock');
            set(gca,'LineWidth',2.0);
            grid on; 
            box off; 
        end
        
 
        function marker = getMarker(obj)
            switch obj.description.type
                case 'pointEstimates'
                    marker = '-';
                case 'CS'
                    marker = ':';
                case 'MC'
                    marker = 'o';
                case 'std'
                    marker = '--';
                otherwise
                    marker = obj.description.type;
            end
                
            
        end
    end
    
end

function [obj,a,b] = genericOperation(obj1,obj2)
            a =  double(obj1);
            b =  double(obj2);
            obj = obj2;
end
