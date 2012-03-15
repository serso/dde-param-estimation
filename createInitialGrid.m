% Creates initial grid: iterating through t from tMin to tMax such that N
% points will be calculated.
%
% @return t - grid for t (array of N elements)
% @return x - grid of values of function f in points of t array (array of N elements)
% @return xWithErrors - grid of values of function f in points of t arrray with normal distribution (mean = 0, sigma = xSigmaErrors) 
% @return minTDist - distance of two closest next points of t grid
% @return maxTDist - distance of two farthest next points of t grid (diameter of grid)

function [t, x, xWithErrors, minTDist, maxTDist] = createInitialGrid (f, N, tMin, tMax, xSigmaError, tSigmaError)

% argument checking

if ( N <= 0)
    throw (MException ('ArgumentCheck:OutOfRange', 'tMax has to be more than tMin'));
end

if ( tMax <= tMin ) 
    throw (MException ('ArgumentCheck:OutOfRange', 'tMax has to be more than tMin'));
end

if ( xSigmaError < 0 ) 
     throw (MException ('ArgumentCheck:OutOfRange', 'xSigmaError has to be not negative'));
end

if ( tSigmaError < 0 ) 
     throw (MException ('ArgumentCheck:OutOfRange', 'tSigmaError has to be not negative'));
end

% function start

t = zeros (N, 1);
x = zeros (N, 1);
xWithErrors = zeros (N, 1);

tStep = (tMax - tMin) / ( N - 1 );

minTDist = Inf;
maxTDist = - Inf;

for i = 1: 1: N
    
    if (i == 1)
        t(i, 1) = tMin;
    elseif (i == N)
        t(N, 1) = tMax;
    else
        t(i, 1) = tMin + (i - 1) * tStep + tStep * randn(1) * tSigmaError;
        
        while ( t(i - 1, 1) > t(i, 1) )
            t(i, 1) = tMin + (i - 1) * tStep + tStep * randn(1) * tSigmaError;
        end
        
    end
    
    if ( i > 1 )
        if ( t(i, 1) - t(i - 1, 1) < minTDist )
            minTDist = t(i, 1) - t(i - 1, 1);
        end
        
        if ( t(i, 1) - t(i - 1, 1) > maxTDist )
            maxTDist = t(i, 1) - t(i - 1, 1);
        end
        
    end
    
    x(i, 1) = f(t(i, 1));
    if ( isnan(x(i, 1)) )
        hrow (MException ('ArgumentCheck:IllegalArgument', 'Delay function returned NaN!'));
    end
    
    xWithErrors(i, 1) =  x(i, 1) + xSigmaError * randn(1) * x(i, 1);
end

end