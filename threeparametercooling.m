function [coefficients, minSSE] = threeparametercooling(x,y)

    [xSorted, sortIndex] = sort(x);
    
    ySorted = y(sortIndex);
    
    minSSE = inf;
    
    n = length(x);
    sumY = sum(y);
    
    xSquared = xSorted.^2;
    xy = xSorted.*ySorted;
    
    for i = 2:n-1
        numGreater = n - i + 1;
        sumYGreater = sum(ySorted(i:end));
        sumXGreater = sum(xSorted(i:end));
        sumXSquaredGreater = sum(xSquared(i:end));
        sumXYGreater = sum(xy(i:end));
        
        b0 = mean(ySorted(1:i-1));
        b1 = (numGreater*sumXYGreater-sumXGreater*sumYGreater)/ ...
                (numGreater*sumXSquaredGreater-sumXGreater*sumXGreater);
        
        N = (n-numGreater)*sumXYGreater*sumXGreater + ...
                sumXSquaredGreater*(numGreater*sumY-n*sumYGreater) + ...
                sumXGreater*sumXGreater*(sumYGreater-sumY);
        D = (n-numGreater)*(numGreater*sumXYGreater-sumXGreater*sumYGreater);
        b2 = N/D;    
        
        
        residuals = y - b0 - b1*(max(x-b2,0));
        
        sse = sum(residuals.^2);
        
        if sse < minSSE
           minSSE = sse;
           coefficients(1) = b0;
           coefficients(2) = b1;
           coefficients(3) = b2;
        end                              
    end                    
end