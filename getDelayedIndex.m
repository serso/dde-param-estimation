function iDelay = getDelayedIndex(t, i, delay)
%GETDELAYEDINDEX finds index of element in t array which value is closest to the
% (t(i) - delay) value. In case if such value cannot be found => returns -1
%
%   iDelay = GETDELAYEDINDEX(t, i, delay):
%   t       vector of availalbe t values
%   i       curent position in t array (t(i) is current position in time),
%           delay index will be found 
%   delay   delay of time which index is needed (t(i) - delay)


% run from the current position i down and find the first element which is
% earlier than current t(i) on more than delay value
j = 0;
while ( i - j >= 1 )
    if ( delay <= t(i) - t(i - j) ) 
        break;
    end
    
    j = j + 1;
end

if ( i - j == 0)
    % unable to find index in t array => we need to use pre history =>
    % return -1
    iDelay = -1;
else 
    iDelay = j;
end

end