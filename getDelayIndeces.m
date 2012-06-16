function result = getDelayIndeces (t, delays, p, options)
%GETDELAYINDECES returns indeces of delayed elements
%
%   result = GETDELAYINDECES(t, delays, p):
%   t       t vector
%   delays  vector of delays, e.g. [0, 1, 14, 25, -1, -2] where negative
%   values are defined according to the next rule: if delays(i) is negative
%   then effective i-th delay is p(-delays(i)), i.e. delays(i) is
%   negative index of delay represented by p vector
%   p   auxiliary vector of delays
%
%   Example: t = [1 2 3 4 5 6 7 8 9], delays = [0 2 7] => result = [1 3 8]

methods = struct('left',1, 'central',2, 'right',3);
method = methods.left;

tMin = t(1);
i = 1;

j0 = 1;

ndelays = length(delays);
result = zeros(ndelays, 1);

%origDelays = delays;

[delays, ascDelays] = getDelays(delays, p);
for delay = delays
    
    if ( options.tEven )
        % t is even => may take 2 elements and count step
        tStep = t(2) - t(1);
        result(i) = floor(delay / tStep) + 1;
        i = i + 1;
    else
        
        if ( delay >= 0 )
            % delay => let's iterate throuh t array and find the index of
            % element which is delayed more than current delay
            
            for j = j0:length(t)
                
                % tj - time offset at the index j
                tj = t(j) - tMin;
                
                if ( tj - delay >= 0 )
                    
                    if ( method == methods.left )
                        if ( j > 1 )
                            result(i) = j - 1;
                        else
                            result(i) = 1;
                        end
                    else
                        % we exceed the delay value => found the index => stop
                        if (j > 1)
                            
                            if ( method == methods.central )
                                % prev_tj - time offset at the index previous to j: j - 1
                                prev_tj = t(j - 1) - tMin;
                                
                                % if previous point is closer to delay use j-1 otherwise j
                                if ( delay - prev_tj < tj - delay )
                                    result(i) = j - 1;
                                else
                                    result(i) = j;
                                end
                            else
                                result(i) = j;
                            end
                        else
                            result(i) = j;
                        end
                    end
                    
                    if ( ascDelays )
                        % if all delays are ascendent sorted we can save last j to start from
                        % it next time
                        j0 = j;
                    else
                        j0 = 1;
                    end
                    
                    i = i + 1;
                    break;
                    
                end
            end
            
        else
            display(delays);
            
            % throw (MException ('AssertionError:ConditionFailed', 'Delay must be positive!'));
            
            % no way
            result(i) = 0;
            i = i + 1;
        end
    end
end

result = result';

end