classdef multivariateTimeSeries < handle
    %multivariateTimeSeries Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private )
        labelsOfTimeSeries = [];  
        tsInColumns = [] ;
        configuration;
    end
    
    methods
        function obj = multivariateTimeSeries
            obj.configuration = configFile;

            obj.tsInColumns  = csvread(obj.configuration.dataFilenameCSV,1);
            obj.tsInColumns  = obj.configuration.prepareRawData(obj.tsInColumns);
            
            obj.labelsOfTimeSeries = readCSVheader(obj.configuration.dataFilenameCSV,obj.countTS);
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
    end
    
end

function firstRow = readCSVheader(dataFilenameCSV,nColumns)
fid       = fopen(dataFilenameCSV);
firstRow = textscan(fid,[repmat('%[^,],',1,nColumns-1) '%[^,\r\n]'], 1);
fclose(fid);
end

