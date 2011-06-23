% Method created x vector with time delays (e.g. x = [ x(t-t1) x(t-t2) x(t - t3) ] where tis are elements of delays vector delays = [ t1, t2, t3 ]  )
%
% @param x - x point at t (NO DELAY)
% @param t - t point
% @param i - current position inj t, x arrays
% @params delays - vector of delays
% @params delayF - handler to function which is responsible for x values
% before t(0)
%
% @return xResult - vector of length(delays) elements: xResult = [ x(t-t1) x(t-t2) x(t - t3) ]
%
function xResult = getDelayedX ( x, t, i, delays, delayF )

xResult = zeros(length(delays), 1);

j = 1;
for delay = delays
    % processing current delay
    delayIndex = getDelayedIndex(t, i, delay);
    
    if ( delayIndex ~= -1 ) 
        % delayed x can be get from from x array
        xResult(j) = x (i - delayIndex);
    else
        % no history is available - use function to get delayed x
        xResult(j) = delayF( -delay + t(i) );
        if ( isnan(xResult(j)) )
            throw (MException ('ArgumentCheck:IllegalArgument', 'Delay function returned NaN!'));
        end
    end
    
    j = j + 1;
end
