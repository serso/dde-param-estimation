function xDelays = getDelayedX ( x, t, i, delays, h )
%GETDELAYEDX return vector x calculated at the points t(i)-τ_1, ... ,
% t(i)-τ_N: [x(t(i)-τ_1), ..., x(t(i)-τ_N)]. Note: τ_i might be equal to 0.
%
%   xDelays = GETDELAYEDX( x, t, i, delays, h ):
%   x       vector of x values
%   t       vector of t values
%   i       curent position in t array (t(i) is current position in time)
%   delays  delays of time which calculation of x is needed
%   h       handler for history function to calculate values of x(t(i)-τ_N)
%   if t(i)-τ_N < t(0) (=min(t))

xDelays = zeros(length(delays), 1);

j = 1;

% delays exceeding current threshold might be calculated only with history
% function
delayThreshold = Inf;
for delay = delays
    % processing current delay
    
    if ( delay <= delayThreshold )
        % delay is less than threshold => try to find index
        delayIndex = getDelayedIndex(t, i, delay);
    else
        % delay is more than threshold => history function must be used
        delayIndex = -1;
    end
    
    if ( delayIndex ~= -1 ) 
        % delayed x can be get from from x array
        xDelays(j) = x (i - delayIndex);
    else
        if ( delayThreshold > delay ) 
            delayThreshold = delay;
        end
        % no history is available - use function to get delayed x
        xDelays(j) = h( -delay + t(i) );
        if ( isnan(xDelays(j)) || isempty(xDelays(j)) )
            display(xDelays(j));
            throw (MException ('ArgumentCheck:IllegalArgument', 'Delay function returned bad value!'));
        end
    end
    
    j = j + 1;
end