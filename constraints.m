% Prepares constraints
%
% @param derivative - handle to derivative function ( =f(x, t, theta) )
% @param x - x point
% @param t - t point
% @param N - dimension of x (without theta)
% @params p - dimension of parameter theta
% @params method - ethod for which constaints are created
%
% @return cResult - vector of inequalities calculated at the point (x, t)
% @return ceqResult - vector of equalitites calculated at the point (x, t)
%
function [cResult, ceqResult, cDerivativeResult, ceqDerivativeResult, cHessResult, ceqHessResult] = constraints (f, grad, hess, x, t, N, p, delays, delayF, method)

theta = x ( N + 1 : N + p );

if ( strcmp(method, 'euler') )
    
    ceq = zeros (N - 1, 1);
    
    ceqElementI = 1;
    if ( ~isempty(grad) )
        ceqElements = zeros(2 * (N - 1) + p * (N - 1), 3);
    else
        ceqDerivative = [];
    end
    
    if ( ~isempty(hess) )
        ceqHess = zeros(N + p, N -1);
    else
        ceqHess = [];
    end
    
    for i = 1:1:(N - 1)
        
        delayedX = getDelayedX(x, t, i, getDelays(delays, theta), delayF);
        
        ceq(i) = - x(i) + x(i + 1) - delta(t, i) * f(delayedX, t(i), theta);
        
        if ( ~isempty(grad) )
            [xDerivative, thetaDerivative] = getDerivatives(grad(delayedX, t(i), theta), 1, p);
            
            for j = 1:1:N
                if ( j == i )
                    [ceqElements, ceqElementI] = addSparseElement(ceqElements, ceqElementI, j, i, (- 1 - delta(t, i) * xDerivative));
                elseif ( j == i + 1)
                    [ceqElements, ceqElementI] = addSparseElement(ceqElements, ceqElementI, j, i, 1);
                end
            end
            
            for j = (N+1):1:(N+p)
                [ceqElements, ceqElementI] = addSparseElement(ceqElements, ceqElementI, j, i, (- delta(t, i) * thetaDerivative(j - N)));
            end
        end
        
        
    end
    
    if ( ~isempty(grad) )
        ceqDerivative = createSparseMatrix(ceqElements);
    end
    
elseif (strcmp(method, 'backward_euler'))
    
    ceq = zeros (N - 1, 1);
    
    ceqElementI = 1;
    if ( ~isempty(grad) )
        ceqElements = zeros(2 * (N - 1) + p * (N - 1), 3);
    else
        ceqDerivative = [];
    end
    
    if ( ~isempty(hess) )
        ceqHess = zeros(N + p, N -1);
    else
        ceqHess = [];
    end
    
    for i = 1:1:(N - 1)
        delayedX = getDelayedX(x, t, i + 1, getDelays(delays, theta), delayF);
        
        ceq(i) = - x(i) + x(i + 1) - delta(t, i) * f ( delayedX, t ( i + 1), theta );
        
        if ( ~isempty(grad) )
            [xDerivative, thetaDerivative] = getDerivatives(grad(delayedX, t(i + 1), theta), 1, p);
            
            for j = 1:1:N
                if ( j == i )
                    [ceqElements, ceqElementI] = addSparseElement(ceqElements, ceqElementI, j, i, -1 );
                elseif ( j == i + 1)
                    [ceqElements, ceqElementI] = addSparseElement(ceqElements, ceqElementI, j, i, (1 - delta(t, i) * xDerivative) );
                end
            end
            
            for j = (N+1):1:(N+p)
                [ceqElements, ceqElementI] = addSparseElement(ceqElements, ceqElementI, j, i, (- delta(t, i) * thetaDerivative(j - N)) );
            end
        end
        
        if ( ~isempty(hess) )
        end
        
    end
    
    if ( ~isempty(grad) )
        ceqDerivative = createSparseMatrix(ceqElements);
    end
    
elseif (strcmp(method, 'box'))
    
    ceq = zeros (N - 1, 1);
    
    ceqElementI = 1;
    if ( ~isempty(grad) )
        ceqElements = zeros(2 * (N - 1) + p * (N - 1), 3);
    else
        ceqDerivative = [];
    end
    
    if ( ~isempty(hess) )
        ceqHess = zeros(N + p, N -1);
    else
        ceqHess = [];
    end
    
    for i = 1:1:(N - 1)
        delayedX_i = getDelayedX(x, t, i, getDelays(delays, theta), delayF);
        delayedX_i_1 = getDelayedX(x, t, i + 1, getDelays(delays, theta), delayF);
        delayedX = (delayedX_i + delayedX_i_1) / 2;
        
        
        deltaT = delta(t, i);
        if ( i > 1 )
            deltaT =  ( deltaT + delta(t, i - 1))/2;
        end
        ti = t(i) + delta(t, i)/2;
        ceq(i) = - x(i) + x(i + 1) - deltaT * f ( delayedX, ti, theta );
        
        if ( ~isempty(grad) )
            [xDerivative, thetaDerivative] = getDerivatives(grad(delayedX, ti, theta), 1, p);
            
            for j = 1:1:N
                if ( j == i )
                    [ceqElements, ceqElementI] = addSparseElement(ceqElements, ceqElementI, j, i, ( - 1 - deltaT * xDerivative) );
                elseif ( j == i + 1)
                    [ceqElements, ceqElementI] = addSparseElement(ceqElements, ceqElementI, j, i, (1 - deltaT * xDerivative) );
                end
            end
            
            for j = (N+1):1:(N+p)
                [ceqElements, ceqElementI] = addSparseElement(ceqElements, ceqElementI, j, i, (- deltaT * thetaDerivative(j - N)) );
            end
        end
        
        if ( ~isempty(hess) )
        end
        
    end
    
    if ( ~isempty(grad) )
        ceqDerivative = createSparseMatrix(ceqElements);
    end
    
elseif (strcmp(method, 'rk4'))
    
    ceq = zeros (N - 1, 1);
    ceqDerivative = [];
    
    for i = 1:1:(N - 1)
        delayedX = getDelayedX(x, t, i, getDelays(delays, theta), delayF);
        
        K1 = f(delayedX, t(i), theta);
        K2 = f(delayedX + delta(t, i) * K1 / 2, t(i) + delta(t, i) / 2, theta);
        K3 = f(delayedX + delta(t, i) * K2 / 2, t(i) + delta(t, i) / 2, theta);
        K4 = f(delayedX + delta(t, i) * K3, t(i) + delta(t, i), theta);
        
        ceq(i) = - x(i) + x(i + 1) - delta(t, i) * ( K1 + 2 * K2 + 2 * K3 + K4 ) / 6 ;
    end
    
else
    throw (MException ('ArgumentCheck:IllegalArgument', 'Unsupported method!'));
end

% no inequlity constraints
cResult = [];
cDerivativeResult = [];
cHessResult =[];

ceqResult = ceq;
ceqDerivativeResult = ceqDerivative;
ceqHessResult = ceqHess;


end