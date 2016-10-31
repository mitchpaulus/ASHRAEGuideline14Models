function [coefficients, minSSE] = threeparametercooling(x,y)

    [xSorted, sortIndex] = sort(x); %Sort the input arrays by increasing x    
    ySorted = y(sortIndex);
    
    minSSE = inf;  %initially set min SSE to arbitrarily high value
    
    %Calculate variable that are unrelated to location of split decision
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
        
        b0 = mean(ySorted(1:i-1));   %EQ. 17
        b1 = (numGreater*sumXYGreater-sumXGreater*sumYGreater)/ ...
                (numGreater*sumXSquaredGreater-sumXGreater*sumXGreater); %EQ. 18
        
        N = (n-numGreater)*sumXYGreater*sumXGreater + ...
                sumXSquaredGreater*(numGreater*sumY-n*sumYGreater) + ...
                sumXGreater*sumXGreater*(sumYGreater-sumY);  %EQ. 20
        D = (n-numGreater)*(numGreater*sumXYGreater-sumXGreater*sumYGreater); %EQ. 21.
        b2 = N/D;  %EQ. 19
                
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
