classdef IRFcollection < IRF 
    %IRFCOLLECTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        function  obj = IRFcollection(IRFmatrix,names,label)
            if nargin ~= 0
                nElements = size(IRFmatrix,1);
                obj(nElements) = IRFcollection;
                for i = 1:nElements
                    currentIRF = IRF(IRFmatrix(i,:),names(i),label);
                    obj(i).Values = currentIRF.Values;
                    obj(i).labelTS = names(i);
                    obj(i).label = label;
                end
            end
        end
        function  obj = setLabel(obj, label)
            nElements = size(obj,2);
            for i = 1 : nElements
                obj(i).label=label;
            end
        end
        function  irfMat =  matrixForm(obj)
            nElements = size(obj,2);
            maxHorizon = size(obj(1).Values,2);
            irfMat = zeros(nElements,maxHorizon);
            for i = 1 : nElements
                irfMat(i,:) = obj(i).Values;
            end
        end
        function  d = double(obj)
            d = obj.matrixForm;
        end
        function  obj = setValues(obj,values)
            nElements = size(obj,2);
            for i = 1 : nElements
                obj(i).Values=values(i,:);
            end
        end
    end
    
end

