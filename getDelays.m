function d = getDelays( delays, theta )
%GETDELAYS calculates effective delay values
%
%   [ci, ce, cij, cej] = GETDELAYS(delays, theta):
%   delays  vector of delays, e.g. [0, 1, 14, 25, -1, -2] where negative
%   values are defined according to the next rule: if delays(i) is negative
%   then effective i-th delay is theta(-delays(i)), i.e. delays(i) is
%   negative index of delay represented by theta vector
%   theta   auxiliary vector of delays

if ( length(delays) < 1 )
    throw (MException ('ArgumentCheck:IllegalArgument', 'Delay vector must contain at least one element'));
end

d = zeros(length(delays), 1);

i = 1;
for delay = delays
    
    if ( delay < 0 )
        d(i) = theta(-delay);
    else
        d(i) = delay;
    end
    
    i = i + 1;
end

d = d';

end
