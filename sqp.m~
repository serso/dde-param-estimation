function [xResult, xStepsResult, timeResult] = sqp (N, lambdaN, fun, grad, hess, c, jc, x0, options, debug, lb, ub)

timerId = tic();

if ( isempty(N) )
    throw (MException ('ArgumentCheck:NotSet', 'N is not set.'));
elseif ( N <= 0 )
    throw (MException ('ArgumentCheck:OutOfRange', 'N is less than 1.'));
end

if ( isempty(grad) )
    throw (MException ('ArgumentCheck:NotSet', 'Gradient must be set.'));
end

if ( isempty(hess) )
    throw (MException ('ArgumentCheck:NotSet', 'Hessian must be set.'));
end

if ( isempty(x0) )
    x0 = zeros(N, 1);
elseif ( length(x0) ~= N )
    throw (MException ('ArgumentCheck:illegalArgument', 'x0 length is not equal to N.'));
end

% if ( isempty(lb) )
%     lb = -Inf * ones ( N, 1 );
% end
% 
% if ( isempty(ub) )
%     ub = Inf * ones ( N, 1 );
% end

maxNumberOfIterations = 5;
tolX = 10^-2;

xResult = x0;

    function [sResult] = direction ( s, x, lambda, hessMatrix )
        sResult = grad(x, lambda)' * s + s' * hessMatrix * s / 2;
    end

    function [cResult, ceqResult] = directionConstraints ( s, x )
        cResult = [];
        if ( ~isempty(c) && ~isempty(jc) )
            ceqResult = c(x) + jc(x) * s;
        else
            ceqResult = [];
        end
    end

for algorithmIteration = 1:2147483647
    
    sLb = lb - xResult;
    sUb = ub - xResult;

    %sLb = -deltaS * ones(N, 1);
    %sUb = deltaS * ones(N, 1);
    
    [~, funGrad] = fun(xResult);
    Ad = jc(xResult);
    lambda = - ( Ad * Ad' ) \ ( Ad * funGrad );
    
    hessMatrix = hess(xResult, lambda);
    
    %% solving the task
    [sResult] = ...
        fmincon( ...
        ...% minimized function
        @(sArg)direction(sArg, xResult, lambda, hessMatrix), ...
        ...% x0
        zeros(N, 1), ...
        ...% A, b (A * x <= b)
        [], [], ...
        ...% Aeq, beq (Aeq * x = beq)
        [], [], ...
        ...% lb, ub (lb <= x < ub)
        sLb, sUb, ...
        @(sArg)directionConstraints(sArg, xResult), ...
        options);
    
    xStepsResult(algorithmIteration, 1:N) = xResult';
    xResult = xResult + sResult;
    
    if ( debug )
        display(sprintf('Iteration: %i', algorithmIteration));
        display(xStepsResult);
    end
    
    if ( algorithmIteration >= maxNumberOfIterations || norm(sResult, inf) < tolX )
        xStepsResult(algorithmIteration + 1, 1:N) = xResult';
        break;
    end
    
end

timeResult = toc(timerId);

if ( debug )
    display(sprintf('Time: %0.5f s', timeResult));
end

end