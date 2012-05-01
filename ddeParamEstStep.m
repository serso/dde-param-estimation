%%
% Function estimate parameters of dde: (*) dx/dt = f (t, x(t), x(t-τ_1), ..., x(t - τ_p)),
% where τ_1, ..., τ_p are constant time delays (which can be also estimated)
% on grid of (t, x) (known points)
%
% @param f - handler to the function f of (*)
% @param fg - handler to the gradient of f (might be [])
% @param fh - handler to the hessian of f (might be [])
% @param t - time grid
% @param x - x grid
% @param N - number of elements in derivative approximation grid
% @param p - number of estimated parameters
% @param delays - vector of delays tau (if delay is estimated - use NUMBER OF ELEMENT in parameter estimation vectors )
% @param h - function for history of x in [tMin - max(tau_i), tMin]
% @param method - derivative approximation method
%
function [xResult, thetaResult, sosResult, exitflagResult, outputResult, lambdaResult, gradResult, hessResult, timeResult] = ...
    ddeParamEstStep ( ...
    f, fg, fh, ...
    t, x, ...
    N, p, ...
    delays, h, ...
    method, ...
    debug, showResult, options, x0, deltaT, thetaLb, thetaUb)

timerId = tic;

obligatoryArgs = 15;

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

%% fill H matrix and g vector for least squares function

% lsf = 
% = ∑ ( y_i - x_i ) ^ 2 =
% = ( y_1 - x_1 ) ^ 2 + ( y_2 - x_2 ) ^ 2 + ... ( y_N - x_iN ) ^ 2 = 
% = ( y_1^2 - 2*y_1*x_1 + x_1^2 ) + ( y_2^2 - 2*y_2*x_2 + x_2^2 ) + ... + ( y_N^2 - 2*y_N*x_N + x_N^2 ) =
% = ∑ y_i^2 - 2 * ∑ y_i*x_i + ∑ x_i^2
% => lsf = x * H * x + g * x + c
% where
% H = one(N, 1)
% g = - 2 [y_1 y_2 ... y_N ]
% c = y^2
%
% array of H matrix elements, where
% first parameter is number of row (i)
% second is number of column (j)
% third is value in i-th row and j-th column
hElements = zeros (N + p, 3);

g = zeros (N + p, 1);

j = 1;
for i = 1 : 1 : ( N + p )
    
    if ( i <= N  )
        
        if ( tUsed(i) > 0 )
            % diagonal elements to 1
            [hElements, ~] = addSparseElement(hElements, i, i, i, 1);
            
            g(i) = - 2 * x(j);
            
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
    display (g);
    display (t);
    display (tGrid);
    display (x);
end

if (strcmp(method, 'rk4') && ~isempty(fg))
    warning ('ArgumentCheck:Warning', 'User-supplied gradient and hessResult are not supported for rk4! Finite differences approximation will be used instead.');
end

if ( ~isempty(fg) && ~strcmp(method, 'rk4') )
    options.optOptions = optimset (options.optOptions, 'GradConstr','on');
else
    options.optOptions = optimset (options.optOptions, 'GradConstr','off');
end

function lResult = hessianF (x_k, lm_k)
    lResult = hessian(x_k, lm_k.eqnonlin, H, fh, tGrid, N, p, method, delays, h);
end

if ( ~isempty(fh) && ~strcmp(method, 'rk4') )
    options.optOptions = optimset (options.optOptions, 'Hessian','user-supplied', 'HessFcn', @hessianF);
else
    options.optOptions = optimset (options.optOptions, 'Hessian','off');
end

options.optOptions = optimset(options.optOptions, 'GradObj', 'on');

%% least squares function

% just to avoid calculations inside fmincon
xConst = x' * x;
H_by_2 = 2 * H;

    function [lsf_, lsfg, lsfh] = lsf ( x )
        lsf_ = x' * H * x + g' * x + xConst;
        lsfg = H_by_2 * x + g;
        lsfh = H_by_2;
    end

    function cej_ = cej(x)
        [~, ~, ~, cej_] = constraints([], fg, [], tGrid, x, N, p, delays, h, method);
    end

    function lsfh = lsf_lh (x, lm_k)
        if ( strcmp(options.hessian_method, 'newton') )
            %lsfh = H_by_2 + hessian(x, lm_k, H, fh, tGrid, N, p, method, delays, h);
            [~, ~, ~, ~, ~, ceh] = constraints([], [], fh, tGrid, x, N, p, delays, h, method, lm_k(N+p:N+N+p-1));
            % todo serso: check the lm_k parameter
            lsfh = H_by_2 + ceh;
        elseif ( strcmp(options.hessian_method, 'gauss-newton') )
            lsfh = H_by_2;
        else
            throw (MException ('ArgumentCheck:IllegalArgument', 'Unsupported method!'));
        end
    end

    function ce_ = ce (x) 
        [~, ce_] = constraints(f, [], [], tGrid, x, N, p, delays, h, method);
    end

    function [ci, ce, cij, cej] = constraints0 (x) 
        [ci, ce, cij, cej] = constraints(f, fg, [], tGrid, x, N, p, delays, h, method);
        cej = cej';
    end


%% solving the task

% sqpTimerId = tic;

if ( options.sqp )
    
    options.sqpOptions.dxmin = options.xTol / 1000;
    
    [solution, lambdaResult, iterations, funCount] = sqp(    ...
        N + p, N - 1, ...
        @(x_k)lsf(x_k), ...
        @(x_k, lm_k)lsf_lh(x_k, lm_k), ...
        @(x_k)ce(x_k), ...
        @(x_k)cej(x_k), ...
        x0, ...
        options.sqpOptions, ...
        debug, ...
        lb, ub);
    
    sosResult = lsf(solution);
    exitflagResult = 1;
    outputResult = struct('iterations', iterations, 'funcCount', funCount);
    gradResult = 1;
    hessResult = 1;
    
else
    
    options.optOptions.TolX = options.xTol / 1000;
    
    [solution, sosResult, exitflagResult, outputResult, lambdaResult, gradResult, hessResult] = ...
        fmincon( ...
        @(x_k)lsf(x_k), ...
        x0, [], [], [], [], ...
        lb, ub, ...
        @(x_k)constraints0(x_k), ...
        options.optOptions);
    
end

% sqpTimeResult = toc(sqpTimerId);

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