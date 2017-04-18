classdef tensorOperations
    %TENSOROPERATIONS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        function convolution = convWithVector(tensor,convDimension,vector)
            shapeVectorInput = size(tensor);
            nDim = size(shapeVectorInput,2);
            
            originalOrder = 1 : nDim;
            permutedOrder = originalOrder;
            permutedOrder(convDimension) = 1;
            permutedOrder(1) = convDimension;
            tensorPermuted = permute(tensor,permutedOrder);
            
            productTensor = tensorOperations.vectorDotTensor(vector, tensorPermuted);
            
            firstColumnIndexShifted = convDimension-1;
            shortOrder = 1 : (nDim-1);
            reverseOrder =  [firstColumnIndexShifted shortOrder(1:(firstColumnIndexShifted-1))   shortOrder(convDimension:end)  ];
            
            convolution = permute(productTensor,reverseOrder);
        end
        
        
        
        function productTensor = vectorDotTensor(vector, tensor)
            shapeVectorInput = size(tensor);
            nRows = shapeVectorInput(1);
            shapeVectorOutput =shapeVectorInput(2:end);
            nCols = prod( shapeVectorOutput);
            
            rowVector = tensorOperations.rowForm(vector);
            matrixForm = reshape(tensor,[nRows,nCols]);
            productVector = rowVector * matrixForm;
            productTensor = reshape(productVector,shapeVectorOutput);
        end
        
        function rowVector = rowForm(vector)
            if isrow(vector)
                rowVector = vector;
            else
                rowVector = vector';
            end
        end
    end
    
end

