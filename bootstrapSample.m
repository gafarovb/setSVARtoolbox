classdef bootstrapSample
    %BOOTSTRAPSAMPLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        values;
    end
    
    methods
        function obj = bootstrapSample(sample)
            obj.values =  sample;
        end
        function stdVector = std(obj)
            stdVector = std(obj.values,1,2) ;
        end
        function nSamples = getN(obj)
            nSamples = size(obj.values, 2) ;
        end
        function sampleStd = standardized(obj)
            sampleStd =  obj.values ./ repmat( obj.std, 1, obj.getN);
        end
        
        function sampleStud = studentized(obj)
            sampleStud = standardized(obj) - repmat(mean(standardized(obj),2),  1, obj.getN);
        end
        function c = quantile(obj,level)
            c = quantile(obj.values,level,2);
        end
    end
    
end

