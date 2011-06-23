function delayIndecesResult = getDelayIndeces (t, delays, theta)

tMin = min(t);
k = 1;
delayIndecesResult = zeros(length(delays), 1);
for delay = getDelays(delays, theta)
    if ( delay > 0 )
        for i = 1:length(t)
            if ( t(i) - tMin >= delay - 0.00000001 )
                delayIndecesResult(k) = i;
                k = k + 1;
                break;
            end
        end
    else
       delayIndecesResult(k) = 0;
       k = k + 1;
    end
end

delayIndecesResult = delayIndecesResult';

end