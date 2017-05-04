classdef IRFcollection < IRF
    %IRFCOLLECTION Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function  obj = IRFcollection(IRFmatrix, TSdescriptions, IRFdescription)

            if nargin ~= 0
                nElements = size(IRFmatrix,1);
                obj(nElements) = IRFcollection;
                for i = 1:nElements
                    obj(i).Values = IRFmatrix(i,:);
                    if ~isempty(TSdescriptions)
                        obj(i).TSdescription = TSdescriptions(:,i);
                    end
                    if ~isempty(IRFdescription)
                        obj(i).description = IRFdescription;
                    end
                    
                end
            end
        end
        function  obj = setValues(obj,values)
            nElements = size(obj,1) * size(obj,2);
            for i = 1 : nElements
                obj(i).Values=values(i,:);
            end
        end
        function  obj = setDescription(obj, description)
            nElements = size(obj,1) * size(obj,2);
            for i = 1 : nElements
                obj(i).description=description;
            end
        end
        
        function  obj = setDescriptionField(obj, fieldName, newString)
            nElements = size(obj,1) * size(obj,2);
            for i = 1 : nElements
                obj(i).description.(fieldName) = newString;
            end
        end
 
        function  irfMat =  matrixForm(obj)
            irfArray =  arrayForm(obj);
            shapeOfArray = size(irfArray);
            nLines = size(irfArray,3);
            irfMat = reshape(irfArray,shapeOfArray(1)*nLines,shapeOfArray(2));
        end
        
        function  irfArray =  arrayForm(obj)
            [nLines, nTS] = size(obj);
            maxHorizon = size(obj(1).Values,2);
            irfArray = zeros(nTS,maxHorizon,nLines);
            for j = 1 :nLines
                for i = 1 : nTS
                    irfArray(i,:,j) = obj(i).Values;
                end
            end
        end
        
        function  d = double(obj)
            d = obj.matrixForm;
        end
        function  disp(obj)
            [nLines, ~] = size(obj);
             for j = nLines:-1:1
                Description(j) = {obj(j,1).description.legend};
            end
             t =  table;
             t.Description = Description';
             disp(t);
        end
         
        function IRFColl =  join(obj,varargin)
            nVarargs = length(varargin);
            switch nVarargs
                case 0
                    IRFColl = obj;
                case 1
                    IRFColl =  [ obj ;  varargin{nVarargs}];
                otherwise
                    IRFColl = [ join(obj,varargin{1:(nVarargs-1)}) ;  varargin{nVarargs}];
            end
        end
                
        function figureArray = plotPanel(obj,printFilePathAndPrefix)
            [nLines, nPlots] = size(obj);
            for j = 1:nPlots
               figureArray(j) = figure;
                for i = 1 : nLines
                    plot(obj(i,j));
                    hold on;
                end
                if nargin>1
                     fn_print(figureArray(j),[printFilePathAndPrefix obj(i,j).getLabelTS]);
                end
            
            end
        end
    end
    
end



function fn_print(h,name)

set(h, 'PaperUnits','centimeters');
set(h, 'Units','centimeters');
pos=get(h,'Position');
set(h, 'PaperSize', [pos(3) pos(4)]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);

print(name,'-dpdf')
print(name,'-depsc2')


end