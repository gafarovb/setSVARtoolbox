classdef multivariateTimeSeries < handle
    %multivariateTimeSeries Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private )
        tsInColumns = [] ;
    end
    properties (Access = public)
        labelsOfTimeSeries = [];
        unitsOfMeasurement = [];
    end
    methods
        function obj = multivariateTimeSeries(tsInColumns, TSdescription )
            obj.tsInColumns  = tsInColumns;
            if nargin<2
                TSdescription = repmat({'UnknownTimeSeries'},1,countTS(obj));
            end
            
            obj.labelsOfTimeSeries  = TSdescription(1,:);
            if (size(TSdescription,1)<2)
                unitsOfMeasurement = repmat({'Units of measurement'},1,countTS(obj));
            else
                unitsOfMeasurement = TSdescription(2,:);
            end
            obj.unitsOfMeasurement = unitsOfMeasurement;
        end
        function n = countTS(obj)
            [~,n] = size( obj.tsInColumns);
        end
        function T = countTimePeriods(obj)
            [T,~] = size( obj.tsInColumns);
        end
        function [Y,X] = getYX(obj,nLags)
            Y = obj.tsInColumns((nLags+1):end,:); %The rows of this matrix are Y_t'
            availableT = size(Y,1);
            laggedY = lagmatrix(obj.tsInColumns,1:1:nLags);
            
            X = [ones(availableT,1), laggedY((nLags+1):end,:)];    %The rows of this matrix are [1,X_{t}']
        end
        function names = getNames(obj)
            names = obj.labelsOfTimeSeries;
        end
        function unitsOfMeasurement = getUnitsOfMeasurement(obj)
            unitsOfMeasurement = obj.unitsOfMeasurement;
        end
        

    end

end



