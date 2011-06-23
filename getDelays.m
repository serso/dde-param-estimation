function delaysResult = getDelays( delays, theta )

delaysResult = zeros(length(delays), 1);

i = 1;
for delay = delays
    if ( delay < 0 )
         delaysResult(i) = theta(-delay);
    else
        delaysResult(i) = delay;
    end
    
    i = i + 1;
end

delaysResult = delaysResult';

end
