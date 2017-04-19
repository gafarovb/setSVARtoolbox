classdef waitBarCustomized < handle
    %WAITBAR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        messageText = 'Please wait';
        maxSteps;
        waitBarHandle;
        startTime;
    end
    
    methods
        function obj = waitBarCustomized(maxSteps)
            obj.maxSteps = maxSteps;
        end
        function setMessage(obj, messageToShow)
            obj.messageText = messageToShow;
        end
        function showProgress(obj, step)
            currentProgress = step/obj.maxSteps;
            switch step
                case 1
                    obj.startTime = tic;
                    obj.waitBarHandle = waitbar(currentProgress , [obj.messageText obj.elapsedTime obj.estimatedTime(step)]  );
                case obj.maxSteps
                    close( obj.waitBarHandle);
                otherwise
                    waitbar( currentProgress, obj.waitBarHandle, [obj.messageText obj.elapsedTime obj.estimatedTime(step)]);
            end
        end
        
        function elapsedTimeString =  elapsedTime( obj)
                    elapsedTimeString = ['; Elapsed time: ' num2str( toc(obj.startTime)) ' sec' ];
        end
        
        function estimatedTimeString = estimatedTime(obj,step)
            elapsedTime = toc(obj.startTime);
            averageTime = elapsedTime/step;
            expTime = averageTime * (obj.maxSteps -step );
            estimatedTimeString = ['; Estimated time: ' num2str(expTime) ' sec' ];
        end
        
    end
    
end
 