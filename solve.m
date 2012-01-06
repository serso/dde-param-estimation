%%
% Function estimate parameters of dde: (*) dx/dt = f (t, x(t), x(t-tau_1), ..., x(t - tau_N)),
% where tau_1, ..., tau_n are constant time delays (which can be also estimated)
% on grid of (t, x)
%
% @param odeF - handler to function f of (*)
% @param odeFGrad - handler to gradient of f odeF (might be [])
% @param odeFHess - handler to hessian of f odeF (might be [])
% @param t - time grid
% @param x - x grid
% @param N - number of elements in derivative approximation grid
% @param p - number of estimated parameters
% @param delays - vector of delays tau (if delay is estimated - use NUMBER OF ELEMENT in parameter estimation vectors )
% @param delayF - function for history of x in [tMin - max(tau_i), tMin]
% @param method - derivative approximation method
%
function [xResult, thetaResult, sosResult, exitflagResult, outputResult, lambdaResult, gradResult, hessResult, timeResult] = ...
    solve ( ...
    odeF, odeFGrad, odeFHess, ...
    t, x, ...
    N, p, ...
    delays, delayF, ...
    method, ...
    debug, showResult, options, x0, thetaLb, thetaUb)

timerId = tic;

obligatoryArgs = 14;

%% checking input arguments

if ( length(t) ~= length(x) )
    throw (MException ('ArgumentCheck:IllegalArgument', 't and x have to have same length'));
end

%%
lb = -Inf * ones(N + p, 1);
if ( nargin > obligatoryArgs )
    lb(N + 1 : N + p ) = thetaLb;
end

ub = Inf * ones(N + p, 1);
if ( nargin > obligatoryArgs )
    ub(N + 1 : N + p ) = thetaUb;
end


%%

if ( debug )
    display('Solving parameter estimation task: ');
    display(method);
end

%% create constraint grid
[tGrid, tUsed] = createGrid(t, N);

%% fill H matrix and f vector for least squares function

% array of H matrix elements, where
% first parameter is number of row (i)
% second is number of column (j)
% third is value in i-th row and j-th column
hElements = zeros (N + p, 3);

f = zeros (N + p, 1);

j = 1;
for i = 1 : 1 : ( N + p )
    
    if ( i <= N  )
        
        if ( tUsed(i) > 0 )
            % diagonal elements to 1
            [hElements, ~] = addSparseElement(hElements, i, i, i, 1);
            
            f(i) = - 2 * x(j);
            
            % we need to iterate through the indeces which were used in
            % tGrid
            for fake = 1 : 1 : tUsed(i)
                j = j + 1;
            end
        else
            [hElements, ~] = addSparseElement(hElements, i, i, i, 0);
        end
        
    else
        [hElements, ~] = addSparseElement(hElements, i, i, i, 0);
    end
    
end

H = createSparseMatrix(hElements);

if ( debug )
    display (H);
    display (f);
    display (t);
    display (tGrid);
    display (x);
end

if (strcmp(method, 'rk4') && ~isempty(odeFGrad))
    warning ('ArgumentCheck:Warning', 'User-supplied gradient and hessResult are not supported for rk4! Finite differences approximation will be used instead.');
end

if ( ~isempty(odeFGrad) && ~strcmp(method, 'rk4') )
    options = optimset (options, 'GradConstr','on');
else
    options = optimset (options, 'GradConstr','off');
end

function lResult = hessianF (xArg, lambdaArg)
    lResult = hessian(xArg, lambdaArg.eqnonlin, H, odeFHess, tGrid, N, p, method, delays, delayF);
end

if ( ~isempty(odeFHess) && ~strcmp(method, 'rk4') )
    options = optimset (options, 'Hessian','user-supplied', 'HessFcn', @hessianF);
else
    options = optimset (options, 'Hessian','off');
end

options = optimset(options, 'GradObj', 'on');

%% least squares function

% just to avoid calculations inside fmincon
xConst = x' * x;

    function [lsfResult, lsfGradResult, lsfHessResult] = lsf ( x )
        H_2 = 2 * H;
        lsfResult = x' * H * x + f' * x + xConst;
        lsfGradResult = (x' * H_2)' + f;
        lsfHessResult = H_2;
    end

    function lsfGradResult = lsfLagrangianGrad (x, lambda)
        Ad =  sqpConstraintsJacobian(tGrid, x, N, p, odeFGrad, method, delays, delayF);
        [~, lsfGradResult] = lsf(x);
        lsfGradResult = lsfGradResult + (lambda' * Ad)';
    end

    function lsfHessResult = lsfLagrangianHess (x, lambda)
        [~, ~, lsfHessResult] = lsf(x);
        lsfHessResult = lsfHessResult + hessian(x, lambda, H, odeFHess, tGrid, N, p, method, delays, delayF);
    end

    function xResult = sqpConstraints (x)
        [~, xResult] = constraints(odeF, [], [], x, tGrid, N, p, delays, delayF, method);
    end

%% solving the task

fminconTimerId = tic;

[solution, sosResult, exitflagResult, outputResult, lambdaResult, gradResult, hessResult] = ...
    fmincon( ...
    @(xArg)lsf(xArg), ...
    x0, [], [], [], [], ...
    lb, ub, ...
    @(xArg)constraints(odeF, odeFGrad, [], xArg, tGrid, N, p, delays, delayF, method), ...
    options);

fminconTimeResult = toc(fminconTimerId);

sqpOptions = optimset('Algorithm', 'interior-point', 'MaxFunEvals', 3000);
sqpOptions = optimset(sqpOptions, 'DerivativeCheck', 'off');
sqpOptions = optimset(sqpOptions, 'FinDiffType', 'central');
sqpOptions = optimset(sqpOptions, 'TolFun', 1e-2);
sqpOptions = optimset(sqpOptions, 'TolCon', 1e-2);
sqpOptions = optimset(sqpOptions, 'TolX', 1e-2);

% sqpTimerId = tic;
% 
% [sqpSolution] = sqp(    ...
%     N + p, N - 1, ...
%     @(xArg)lsf(xArg), ...
%     @(xArg, lambda)lsfLagrangianGrad(xArg, lambda), ...
%     @(xArg, lambda)lsfLagrangianHess(xArg, lambda), ...
%     @(xArg)sqpConstraints(xArg), ...
%     @(xArg)sqpConstraintsJacobian(tGrid, xArg, N, p, odeFGrad, method, delays, delayF), ...
%     x0, ...
%     sqpOptions, ...
%     debug, ...
%     lb, ub);
% 
% sqpTimeResult = toc(sqpTimerId);
% 
% thetaResult = solution( N + 1 : N + p );
% sqpThetaResult = sqpSolution( N + 1 : N + p );
% 
% display('ThetaResult - SqpThetaResult');
% display(thetaResult - sqpThetaResult);
% 
% display(sprintf('Time for fmincon: %0.3f s', fminconTimeResult));
% display(sprintf('Time for sqp: %0.3f s', sqpTimeResult));

%% preparing result
thetaResult = solution( N + 1 : N + p );
% sqpThetaResult = sqpSolution( N + 1 : N + p );

% display('ThetaResult - SqpThetaResult');
% display(thetaResult - sqpThetaResult);

if ( showResult )
    figure;
    grid on;
    hold on;
    title(sprintf('Method: %s\nApproximation grid: %i\nTime: %0.3f s', method, N, toc(timerId)));
    plot (t, x, '-b');
    plot (interpolate(t, N, 'spline'), solution(1:N), 'or');
end

xResult = zeros(length(t), 1);
j = 1;
for i = 1 : 1 : N
    if ( tUsed(i) > 0 )
        
        xResult(j) = solution(i);
        
        % we need to iterate through the indeces which were used in
        % tGrid
        for fake = 1 : 1 : tUsed(i)
            j = j + 1;
        end
    end
end

if ( debug )
    display(xResult);
    display(sosResult);
end

timeResult = toc(timerId);

if (debug || showResult)
    display(sprintf('Method: %s', method));
    display(sprintf('Approximation grid: %i', N));
    display(sprintf('Time: %0.3f s', timeResult));
    
    display(outputResult);
    display(thetaResult);
end

end