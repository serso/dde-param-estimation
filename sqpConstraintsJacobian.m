function [jacobianResult] = sqpConstraintsJacobian (t, x, N, p, odeFGrad, method, delays, delayF)

theta = x ( N + 1 : N + p );

delayIndeces = getDelayIndeces(t, delays, theta);

delayElements = 0;
for delayIndex = delayIndeces
    if (delayIndex > 0)
        delayElements = delayElements + (N - delayIndex - 1);
    end
end

jElements = zeros(2 * (N - 1) + p * (N - 1) + delayElements, 3);


if ( strcmp(method, 'euler') )
    throw (MException ('ArgumentCheck:IllegalArgument', 'Unsupported method!'));
    
elseif ( strcmp(method, 'backward_euler') )
    prevC_i = [];
    j = 1;
    
    for i = 1: 1 :(N - 1)
        delayedX = getDelayedX(x, t, i + 1, getDelays(delays, theta), delayF);
        
        [xDerivative, thetaDerivative] = getDerivatives(odeFGrad(delayedX, t(i + 1), theta), 1, p);
        
        val = delta(t, i) * xDerivative;
        d_i = -1;
        c_i = 1 - val;
        e_i = - delta(t, i) * thetaDerivative;
        delay_i = - val;
        
        [jElements, j] = addSparseElement(jElements, j, i, i, d_i);
        if ( ~isempty(prevC_i) )
            [jElements, j] = addSparseElement(jElements, j, i-1, i, prevC_i);
        end
        
        for delayIndex = delayIndeces
            if delayIndex > 0
                if ( i - delayIndex > 0 )
                   [jElements, j] = addSparseElement(jElements, j, i, i - delayIndex, delay_i);
                end
            end
        end
        
        prevC_i = c_i;
        
        for eIndex = 1:p
            [jElements, j] = addSparseElement(jElements, j, i, N + eIndex, e_i(eIndex));
        end
    end
    
    if ( ~isempty(prevC_i) )
        [jElements, j] = addSparseElement(jElements, j, N-1, N, prevC_i);
    end
    
elseif (strcmp(method, 'rk4'))
    %             not suppported
    throw (MException ('ArgumentCheck:IllegalArgument', 'Unsupported method!'));
else
    throw (MException ('ArgumentCheck:IllegalArgument', 'Unsupported method!'));
end

jacobianResult = createSparseMatrix(jElements);
% display(full(jacobianResult));

end