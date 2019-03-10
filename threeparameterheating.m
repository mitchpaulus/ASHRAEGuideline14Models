function [coefficients, minSSE] = threeparameterheating(x,y)
% THREEPARAMETERHEATING
%   USAGE: 
%     [coefficients, minSSE] = threeparameterheating(x, y)
%
%   x: 1-D vector of x values
%   y: 1-D vector of y values
%
%   coefficients: vector consisting of [constant, slope, changepoint]
%   minSSE: scalar minimum of sum of squared errors
%   
%   Notes: This function employs the algorithm described in 
%
%   Mitchell T. Paulus. Algorithm for explicit solution to the three parameter
%   linear change point regression model. Science and Technology for the Built
%   Environment, 2016
%
%   You can find a copy at:
%   https://www.researchgate.net/publication/312051080_Algorithm_for_explicit_solution_to_the_three_parameter_linear_change_point_regression_model

    [xSorted, sortIndex] = sort(x); %Sort the input arrays by increasing x    
    ySorted = y(sortIndex);
    
    minSSE = inf;  %initially set min SSE to arbitrarily high value
    
    %Calculate variable that are unrelated to location of split decision
    n = length(x);
    sumY = sum(y);    
    xSquared = xSorted.^2;
    xy = xSorted.*ySorted;
    
    for m = 3:n
        L = m - 1; %Using capital L because lowercase l looks too close to 1. 
        sumYLower = sum(ySorted(1:L));  
        sumXLower = sum(xSorted(1:L));
        sumXSquaredLower = sum(xSquared(1:L));
        sumXYLower = sum(xy(1:L));
        
        b0 = mean(ySorted(m:n));   %EQ. 27
        b1 = (sumXLower * sumYLower - L * sumXYLower) / ...
                (L * sumXSquaredLower - sumXLower * sumXLower); %EQ. 28
        
        N = (n - L) * (sumXLower * sumXYLower) + L * sumXSquaredLower * sumY - ...
                sumY * sumXLower * sumXLower - n * sumXSquaredLower * sumYLower + ...
                sumYLower * sumXLower * sumXLower;  %EQ. 20
        D = (n - L) * (L * sumXYLower - sumXLower  * sumXLower); %EQ. 31.
        b2 = N/D;  %EQ. 29
                
        residuals = y - b0 - b1*(max(b2 - x,0));
        
        sse = sum(residuals.^2);
        
        if sse < minSSE
           minSSE = sse;
           coefficients(1) = b0;
           coefficients(2) = b1;
           coefficients(3) = b2;
        end        
        
    end                    
end
