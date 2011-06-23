% hessian of lagrangian for fmincon

%lambda.eqnonlin
function [lResult] = hessian ( x, lambda, H, odeFHess, t, N, p, method, delays, delayF )

lResult = 2 * H;

theta = x ( N + 1 : N + p );

if ( strcmp(method, 'euler') )
    
    for k = 1:1:(N-1)
        out = deal(odeFHess(getDelayedX(x, t, k, getDelays(delays, theta), delayF), t(k), theta));
        lResult(k, k) = lResult(k, k) - lambda(k) * delta(t, k ) * out(1, 1);
    end
    
elseif ( strcmp(method, 'backward_euler') )
    
    for k = 2:1:N
        out = deal(odeFHess(getDelayedX(x, t, k, getDelays(delays, theta), delayF), t(k), theta));
        lResult(k, k) = lResult(k, k) - lambda(k - 1) * delta(t, k - 1 ) * out(1, 1);
    end
    
elseif (strcmp(method, 'box'))
    
    for k = 1:1:N
        if ( k > 1 )
            delayedX_1 = getDelayedX(x, t, k - 1, getDelays(delays, theta), delayF);
            delayedX_2 = getDelayedX(x, t, k, getDelays(delays, theta), delayF);
            delayedX = (delayedX_1 + delayedX_2) / 2;
            out = deal(odeFHess(delayedX, t(k - 1) + delta(t, k - 1 ) / 2 , theta));
            
            deltaT = delta(t, k - 1 );
            if (k > 2)
                deltaT = ( deltaT + delta(t, k - 2 )) / 2;
            end
            
            val = - lambda(k - 1) * deltaT * out(1, 1);
            lResult(k, k) = lResult(k, k) + val;
            lResult(k, k - 1) = lResult(k, k - 1) + val;
            lResult(k - 1, k) = lResult(k - 1, k) + val;
        end
        
        if ( k < N )
            delayedX_1 = getDelayedX(x, t, k + 1, getDelays(delays, theta), delayF);
            delayedX_2 = getDelayedX(x, t, k, getDelays(delays, theta), delayF);
            delayedX = (delayedX_1 + delayedX_2) / 2;
            out = deal(odeFHess(delayedX, t(k) + delta(t, k ) / 2 , theta));
            
            deltaT = delta(t, k );
            if (k < N - 1)
                deltaT = ( deltaT + delta(t, k + 1 )) / 2;
            end
            val = - lambda(k) * deltaT * out(1, 1);
            lResult(k, k) = lResult(k, k) + val;
            lResult(k, k + 1) = lResult(k, k + 1) + val;
            lResult(k + 1, k) = lResult(k + 1, k) + val;
        end
    end
    
elseif (strcmp(method, 'rk4'))
    %             not suppported
else
    throw (MException ('ArgumentCheck:IllegalArgument', 'Unsupported method!'));
end

for k = (N+1):1:(N+p)
    for l = (N+1):1:(N+p)
        kl = 0;
        for index = 1:1:(N-1)
            out = deal(odeFHess(getDelayedX(x, t, index, getDelays(delays, theta), delayF), t(index), theta));
            kl = kl - lambda(index) * delta(t, index ) * out(k-N, l-N);
        end
        lResult(k, l) = lResult(k, l) + kl;
    end
end

end