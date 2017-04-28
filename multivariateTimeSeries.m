classdef multivariateTimeSeries < handle
    %multivariateTimeSeries Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private )
        labelsOfTimeSeries = [];
        tsInColumns = [] ;
    end
    
    methods
        function obj = multivariateTimeSeries(tsInColumns, labelsOfTimeSeries )
            obj.tsInColumns  = tsInColumns;
            obj.labelsOfTimeSeries  = labelsOfTimeSeries;
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

    end

end



