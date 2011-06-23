% Method returns index of delaeyd element
%
% @param t - t array
% @param i - current position inj t arrays
% @params delay - current delay
% @params delayF - handler to function which is responsible for x values
% before t(0)
%
% @return xResult - vector of length(delays) elements: xResult = [ x(t-t1) x(t-t2) x(t - t3) ]
%
function iResult = getDelayedIndex(t, i, delay)

j = 0;
while ( i - j >= 1 && ( delay - t(i) + t(i - j) > 0.00000001 ) )
    j = j + 1;
end

if ( i - j == 0)
    iResult = -1;
else 
    iResult = j;
end

end