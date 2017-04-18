classdef bootstrapSample
    %BOOTSTRAPSAMPLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        values;
        sampleDimension;
    end
    
    methods
        function obj = bootstrapSample(sample, sampleDimension)
            obj.values =  sample;
            if nargin>1
                obj.sampleDimension = sampleDimension;
            else
                obj.sampleDimension = 2;
            end
        end
        function stdVector = std(obj)
            divideByN = 1;
            
            stdVector = std(obj.values, divideByN, obj.sampleDimension) ;
        end
        function nSamples = getN(obj)
            nSamples = size(obj.values, obj.sampleDimension) ;
        end
        function sampleStd = standardized(obj)
            sampleStd =  obj.values ./ repmat( obj.std, 1, obj.getN);
        end
        
        function sampleStud = studentized(obj)
            sampleStud = standardized(obj) - repmat(mean(standardized(obj), obj.sampleDimension),  1, obj.getN);
        end
        function c = quantile(obj,level)
            c = quantile(obj.values,level, obj.sampleDimension);
        end
    end
    
end

