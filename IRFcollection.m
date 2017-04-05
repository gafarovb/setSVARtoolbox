classdef IRFcollection < IRF
    %IRFCOLLECTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        function  obj = IRFcollection(IRFmatrix,names,label)
            if nargin ~= 0
                n = size(IRFmatrix,1);
                obj(n) = IRFcollection;
                for iVar = 1:n
                    currentIRF = IRF(IRFmatrix(iVar,:),names(iVar),label);
                    obj(iVar).Values = currentIRF.Values;
                    obj(iVar).labelTS = names(iVar);
                    obj(iVar).label = label;
                end
            end
        end
    end
    
end

