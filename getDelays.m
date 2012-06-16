function [delaysResult, ascDelays] = getDelays( delays, p )
%GETDELAYS calculates effective delay values
%
%   [delaysResult, ascDelays] = GETDELAYS(delays, p):
%   delays  vector of delays, e.g. [0, 1, 14, 25, -1, -2] where negative
%   values are defined according to the next rule: if delays(i) is negative
%   then effective i-th delay is p(-delays(i)), i.e. delays(i) is
%   negative index of delay represented by p vector
%   p   auxiliary vector of delays

if ( length(delays) < 1 )
    throw (MException ('ArgumentCheck:IllegalArgument', 'Delay vector must contain at least one element'));
end

ascDelays = true;

delaysResult = zeros(length(delays), 1);

prevDelay = -inf;
i = 1;
for delay = delays
    
    if ( delay < 0 )
        delaysResult(i) = p(-delay);
    else
        delaysResult(i) = delay;
    end
    
    if ( prevDelay > delaysResult(i) )
        ascDelays = false;
    end
    
    % update step
    prevDelay = delaysResult(i);
    i = i + 1;
end

delaysResult = delaysResult';

end
