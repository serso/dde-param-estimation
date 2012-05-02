function [ci, ce, cij, cej, cih, ceh] = constraints ( ... 
    f, fg, fh,      ...
    t, x,           ...
    N, p,           ...
    delays, h,      ... 
    method, lm)
%CONSTRAINTS calculates constraints in point x for the DDE parameters
%estimation problem: (*) dx/dt = f (t, x(t), x(t-τ_1), ..., x(t - τ_N))
%
%   [ci, ce, cij, cej] = CONSTRAINTS(f, fg, hess, x, t, N, p, delays, h, method) calculates next values:
%   ci  vector of inequalities calculated in the point (t, x)
%   ce  vector of equalitites calculated in the point (t, x)
%   cij jacobian matrix of inequalities in the point (t, x)
%   cej jacobian matrix of equalities in the point (t, x)
%   cih hessian matrix of inequalities in the point (t, x)
%   ceh hessian matrix of equalities in the point (t, x)
%
%
% @param f - handler to the function f of (*) (might be [])
% @param fg - handler to the gradient of f (might be [])
% @param fh - handler to the hessian of f (might be [])
% @param t - time grid
% @param x - x grid
% @param N - number of elements in integration grid
% @param p - number of estimated parameters
% @param delays - vector of delays tau (if delay is estimated - use NUMBER OF ELEMENT in parameter estimation vectors )
% @param h - function for history of x in [tMin - max(τ_i), tMin]
% @param method - derivative approximation method
%

%% prepare initial parameters

% extract theta
theta = x ( N + 1 : N + p );

% if no f handler is provided => no need for calculations of ce => set empty
if (~isempty(f))
    % allocate ce array
    ce = zeros (N - 1, 1);
else
    ce = [];
end

% prevailing index for cej
cej_i = 1;


% if no fg handler is provided => no need for calculations of cej => set
% empty
if ( ~isempty(fg) )
    % count total value of non zero element to effectively allocate memory
    delayIndeces = getDelayIndeces(t, delays, theta);
    
    delayElements = 0;
    for delayIndex = delayIndeces
        if (delayIndex > 1)
            % NOTE: diagonal is already counted in the code below
            
            % each delay has a sub diagonal => number of not zero elements
            % can be retrieved by the next formula
            if ( strcmp(method, 'euler') )
                delayElements = delayElements + (N - delayIndex - 1);
            elseif (strcmp(method, 'backward_euler'))
                if ( delayIndex == 2 ) 
                    display(delayIndex)
                    throw (MException ('ArgumentCheck:IllegalArgument', 'Check!'));
                end
                
                delayElements = delayElements + (N - (delayIndex + 1) - 1);
            elseif (strcmp(method, 'box'))
                delayElements = delayElements + (N - delayIndex - 1);
            else
                throw (MException ('ArgumentCheck:IllegalArgument', 'Unsupported method!'));
            end            
        end
    end
    
    % allocating memory
    cejArray = zeros(2 * (N - 1) + p * (N - 1) + delayElements, 3);
else
    cej = [];
end

% prevailing index for ceh
ceh_i = 1;

if ( ~isempty(fh) )
    % count total value of non zero element to effectively allocate memory
    delayIndeces = getDelayIndeces(t, delays, theta);
    
    delayElements = 0;
    for delayIndex = delayIndeces
        if (delayIndex > 1)
            % NOTE: diagonal is already counted in the code below
            
            % each delay has a sub diagonal => number of not zero elements
            % can be retrieved by the next formula
            if ( strcmp(method, 'euler') )
                throw (MException ('ArgumentCheck:IllegalArgument', 'Unsupported method!'));
            elseif (strcmp(method, 'backward_euler'))
                throw (MException ('ArgumentCheck:IllegalArgument', 'Unsupported method!'));
            elseif (strcmp(method, 'box'))
                throw (MException ('ArgumentCheck:IllegalArgument', 'Unsupported method!'));
            else
                throw (MException ('ArgumentCheck:IllegalArgument', 'Unsupported method!'));
            end            
        end
    end
    
    % allocating memory
    cehArray = zeros((N - 1) + 2 * p * (N - 1) + p * p + delayElements, 3);
    
    if ( ~isempty(lm) )
        lm = ones(N-1, 1);
    end
else
   ceh = [];
end

effectiveDelays = getDelays(delays, theta);

%% creates output matrices

if ( strcmp(method, 'euler') )
    
    for i = 1:1:(N - 1)
        
        delayedX = getDelayedX(x, t, i, effectiveDelays, h);
        
        if ( ~isempty(f) )
            ce(i) = - x(i) + x(i + 1) - delta(t, i) * f(delayedX, t(i), theta);
        end
        
        if ( ~isempty(fg) )
            
            [xd, thetad] = getDerivatives(fg(delayedX, t(i), theta), 1, p);
            
            [cejArray, cej_i] = addSparseElement(cejArray, cej_i, i, i, (- 1 - delta(t, i) * xd));
            [cejArray, cej_i] = addSparseElement(cejArray, cej_i, i, i+1, 1);
            
            for j = (N+1):1:(N+p)
                [cejArray, cej_i] = addSparseElement(cejArray, cej_i, i, j, (- delta(t, i) * thetad(j - N)));
            end
            
            for delayIndex = delayIndeces
                if delayIndex > 1
                    % todo serso: check if col is correct
                    col = i - delayIndex;
                    if ( col > 0 )
                        [cejArray, cej_i] = addSparseElement(cejArray, cej_i, i, col, - delta(t, i) * xd);
                    end
                end
            end
        end
        
        if (~isempty(fh)) 
            % fhi - function hessian at the point (t_i, x_i)
            fhi = deal(fh(getDelayedX(x, t, i, getDelays(delays, theta), h), t(i), theta));
            
            [cehArray, ceh_i] = addSparseElement(cehArray, ceh_i, i, i, - lm(i) * delta(t, i) * fhi(1, 1)); 
            for j = 1:p
                [cehArray, ceh_i] = addSparseElement(cehArray, ceh_i, i, N + j, - lm(i) * delta(t, i) * fhi(1, 1 + j));
                [cehArray, ceh_i] = addSparseElement(cehArray, ceh_i, N + j, i, - lm(i) * delta(t, i) * fhi(1 + j, 1));
                [cehArray, ceh_i] = addSparseElement(cehArray, ceh_i, N + j, N + j, - lm(i) * delta(t, i) * fhi(1 + j, 1 + j));
            end
        end
        
    end
    
    if ( ~isempty(fg) )
        cej = createSparseMatrix(cejArray);
    end
    
    if ( ~isempty(fh) )
        ceh = createSparseMatrix(cehArray);
    end
    
elseif (strcmp(method, 'backward_euler'))
    
    for i = 1:1:(N - 1)
        
        delayedX = getDelayedX(x, t, i + 1, effectiveDelays, h);
        
        if ( ~isempty(f) )
            ce(i) = - x(i) + x(i + 1) - delta(t, i) * f ( delayedX, t ( i + 1), theta );
        end;
        
        if ( ~isempty(fg) )
            [xd, thetad] = getDerivatives(fg(delayedX, t(i + 1), theta), 1, p);
            
            [cejArray, cej_i] = addSparseElement(cejArray, cej_i, i, i, -1 );
            [cejArray, cej_i] = addSparseElement(cejArray, cej_i, i, i+1, (1 - delta(t, i) * xd) );
            
            for j = (N+1):1:(N+p)
                [cejArray, cej_i] = addSparseElement(cejArray, cej_i, i, j, (- delta(t, i) * thetad(j - N)) );
            end
            
            for delayIndex = delayIndeces
                if delayIndex > 1
                    % todo serso: check if col is correct (another variants: col = i - delayIndex + 1)
                    col = i - delayIndex + 1;
                    if ( col > 0 )
                        [cejArray, cej_i] = addSparseElement(cejArray, cej_i, i, col, - delta(t, i) * xd);
                    end
                end
            end
        end
        
        if (~isempty(fh)) 
            % fhi - function hessian at the point (t_i, x_i)
            fhi = deal(fh(getDelayedX(x, t, i, getDelays(delays, theta), h), t(i), theta));
            
            tmp = - lm(i) * delta(t, i) * fhi(1, 1);
            [cehArray, ceh_i] = addSparseElement(cehArray, ceh_i, i+1, i, tmp); 
            [cehArray, ceh_i] = addSparseElement(cehArray, ceh_i, i, i + 1, tmp); 
            for j = 1:p
                [cehArray, ceh_i] = addSparseElement(cehArray, ceh_i, i, N + j, - lm(i) * delta(t, i) * fhi(1, 1 + j));
                [cehArray, ceh_i] = addSparseElement(cehArray, ceh_i, N + j, i, - lm(i) * delta(t, i) * fhi(1 + j, 1));
                [cehArray, ceh_i] = addSparseElement(cehArray, ceh_i, N + j, N + j, - lm(i) * delta(t, i) * fhi(1 + j, 1 + j));
            end
        end
        
    end
    
    if ( ~isempty(fg) )
        cej = createSparseMatrix(cejArray);
    end
    
    if ( ~isempty(fh) )
        ceh = createSparseMatrix(cehArray);
    end
    
elseif (strcmp(method, 'box'))
    
    for i = 1:1:(N - 1)
        delayedX_i = getDelayedX(x, t, i, effectiveDelays, h);
        delayedX_i_1 = getDelayedX(x, t, i + 1, effectiveDelays, h);
        delayedX = (delayedX_i + delayedX_i_1) / 2;
        
        
%         if ( i > 1 )
%             deltaT =  ( deltaT + delta(t, i - 1))/2;
%         else
            deltaT = delta(t, i);
%         end
        ti = t(i) + (t(i+1) - t(i)) / 2;
        
        if ( ~isempty(f) )
            ce(i) = - x(i) + x(i + 1) - deltaT * f ( delayedX, ti, theta );
        end
        
        if ( ~isempty(fg) )
            [xd, thetad] = getDerivatives(fg(delayedX, ti, theta), 1, p);
            
            [cejArray, cej_i] = addSparseElement(cejArray, cej_i, i, i, ( - 1 - deltaT * xd / 2) );
            [cejArray, cej_i] = addSparseElement(cejArray, cej_i, i, i+1, (1 - deltaT * xd / 2) );
            
            for j = (N+1):1:(N+p)
                [cejArray, cej_i] = addSparseElement(cejArray, cej_i, i, j, (- deltaT * thetad(j - N)) );
            end
            
            for delayIndex = delayIndeces
                if delayIndex > 1
                    % todo serso: check if col is correct (another variants: col = i - delayIndex + 1)
                    col = i - delayIndex;
                    if ( col > 0 )
                        [cejArray, cej_i] = addSparseElement(cejArray, cej_i, i, col, - delta(t, i) * xd / 2);
                    end
                end
            end
        end
        
    end
    
    if ( ~isempty(fg) )
        cej = createSparseMatrix(cejArray);
    end
    
elseif (strcmp(method, 'rk4'))
    
    cej = [];
    
    for i = 1:1:(N - 1)
        delayedX = getDelayedX(x, t, i, effectiveDelays, h);
        
        if ( ~isempty(f) )
            K1 = f(delayedX, t(i), theta);
            K2 = f(delayedX + delta(t, i) * K1 / 2, t(i) + delta(t, i) / 2, theta);
            K3 = f(delayedX + delta(t, i) * K2 / 2, t(i) + delta(t, i) / 2, theta);
            K4 = f(delayedX + delta(t, i) * K3, t(i) + delta(t, i), theta);
            
            ce(i) = - x(i) + x(i + 1) - delta(t, i) * ( K1 + 2 * K2 + 2 * K3 + K4 ) / 6 ;
        end;
    end
    
else
    throw (MException ('ArgumentCheck:IllegalArgument', 'Unsupported method!'));
end

% no inequlity constraints
ci = [];
cij = [];
cih = [];

if ( ~isempty(cej) )
  %  spy(cej);
end

end