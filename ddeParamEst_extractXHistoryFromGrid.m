function [tResult, xResult, minTDistResult, xHistoryResult] = ddeParamEst_extractXHistoryFromGrid(t, x, xHistoryH, maxDelayParam)

maxDelay = -1;
if ( ~isempty(maxDelayParam) )
    if ( maxDelayParam > maxDelay )
        maxDelay = maxDelayParam;
    end
end

if ( isempty(xHistoryH) && maxDelay > 0 )
    i = 1;
    while ( t(i) - t (1) < maxDelay )
        i = i + 1;
    end
    
    tForHistory = t(1:i+1, 1);
    xForHistory = x(1:i+1, 1);
    
    %xHistoryResult = @(t) interp1(tForHistory, xForHistory, t, 'spline', 'extrap');
    
    % using quick interpolation (linear) is more efficient
    xHistoryResult = @(t) interp1q(tForHistory, xForHistory, t);

    tResult = t(i+1:end, 1);
    xResult = x(i+1:end, 1);
else
    xHistoryResult = xHistoryH;
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