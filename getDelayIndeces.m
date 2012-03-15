function result = getDelayIndeces (t, delays, theta)
%GETDELAYINDECES returns indeces of delayed elements
%
%   result = GETDELAYINDECES(t, delays, theta):
%   t       t vector
%   delays  vector of delays, e.g. [0, 1, 14, 25, -1, -2] where negative
%   values are defined according to the next rule: if delays(i) is negative
%   then effective i-th delay is theta(-delays(i)), i.e. delays(i) is
%   negative index of delay represented by theta vector
%   theta   auxiliary vector of delays
%
%   Example: t = [1 2 3 4 5 6 7 8 9], delays = [0 2 7] => result = [1 3 8]

tMin = t(1);
i = 1;

result = zeros(length(delays), 1);

effectiveDelays = getDelays(delays, theta);
for delay = effectiveDelays 
    if ( delay >= 0 )
        % delay => let's iterate throuh t array and find the index of
        % element which is delayed more than current delay
        
        for j = 1:length(t)
            if ( t(j) - tMin >= delay )
                
                % we exceed the delay value => found the index => stop
                result(i) = j;
                i = i + 1;
                break;
                
            end
        end
        
    else
       display(effectiveDelays);
       
       throw (MException ('AssertionError:ConditionFailed', 'Delay must be positive!'));
       
       % no way
       % result(i) = 0;
       % i = i + 1;
    end
end

result = result';

end