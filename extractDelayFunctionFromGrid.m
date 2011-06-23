function [tResult, xResult, minTDistResult, delayFResult] = extractDelayFunctionFromGrid(t, x, delayF, maxDelayParam)

maxDelay = -1;
if ( ~isempty(maxDelayParam) )
    if ( maxDelayParam > maxDelay )
        maxDelay = maxDelayParam;
    end
end

if ( isempty(delayF) && maxDelay > 0 )
    i = 1;
    while ( t(i) - t (1) < maxDelay )
        i = i + 1;
    end
    
    tForDelay = t(1:(i+1), 1);
    xForDelay = x(1:(i+1), 1);
    
    %delayFResult = @(t) interp1(tForDelay, xForDelay, t, 'spline', 'extrap');
    
    % using quick interpolation (linear) is more efficient
    delayFResult = @(t) interp1q(tForDelay, xForDelay, t);

    tResult = t((i+1):length(t), 1);
    xResult = x((i+1):length(x), 1);
else
    delayFResult = delayF;
    tResult = t;
    xResult = x;
end

minTDistResult = Inf;
for i = 1:length(t)
    if ( i > 1 )
        minTDistResult = min ([minTDistResult t(i)-t(i-1)]);
    end
end


end