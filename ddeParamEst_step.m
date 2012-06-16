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
function [xResult, pResult, sosResult, exitflagResult, infoResult, lmResult] = ...
    ddeParamEst_step ( ...
    f, fg, fh, ...
    t, x, ...
    n, np, ...
    delays, xHistory, ...
    method, ...
    debug, showResult, options, x0, deltaT, pLb, pUb)

obligatoryArgs = 15;

%% checking input arguments

if ( length(t) ~= length(x) )
    throw (MException ('ArgumentCheck:IllegalArgument', 't and x have to have same length'));
end

%%
lb = -Inf * ones(n + np, 1);
if ( nargin > obligatoryArgs )
    lb(n + 1 : n + np ) = pLb;
end

ub = Inf * ones(n + np, 1);
if ( nargin > obligatoryArgs )
    ub(n + 1 : n + np ) = pUb;
end


%%

if ( debug )
    display('Solving parameter estimation task: ');
    display(method);
end

%% create constraint grid
[tGrid, tUsed, options.tEven] = ddeParamEst_createGrid(t, n);

%% fill H matrix and g vector for least squares function

% lsf = 
% = ∑ ( y_i - x_i ) ^ 2 =
% = ( y_1 - x_1 ) ^ 2 + ( y_2 - x_2 ) ^ 2 + ... ( y_N - x_iN ) ^ 2 = 
% = ( y_1^2 - 2*y_1*x_1 + x_1^2 ) + ( y_2^2 - 2*y_2*x_2 + x_2^2 ) + ... + ( y_N^2 - 2*y_N*x_N + x_N^2 ) =
% = ∑ y_i^2 - 2 * ∑ y_i*x_i + ∑ x_i^2
% => lsf = x * H * x + g * x + c
% where
% H = one(n, 1)
% g = - 2 [y_1 y_2 ... y_N ]
% c = y^2
%
% array of H matrix elements, where
% first parameter is number of row (i)
% second is number of column (j)
% third is value in i-th row and j-th column
hElements = zeros (n + np, 3);

g = zeros (n + np, 1);

j = 1;
for i = 1 : 1 : ( n + np )
    
    if ( i <= n  )
        
        if ( tUsed(i) > 0 )
            % diagonal elements to 1
            [hElements, ~] = utils.addSparseElement(hElements, i, i, i, 1);
            
            g(i) = - 2 * x(j);
            
            % we need to iterate through the indeces which were used in
            % tGrid
            for fake = 1 : 1 : tUsed(i)
                j = j + 1;
            end
        else
            [hElements, ~] = utils.addSparseElement(hElements, i, i, i, 0);
        end
        
    else
        [hElements, ~] = utils.addSparseElement(hElements, i, i, i, 0);
    end
    
end

H = utils.createSparseMatrix(hElements);

% if ( debug )
%     display (H);
%     display (g);
%     display (t);
%     display (tGrid);
%     display (x);
% end

if (strcmp(method, 'rk4') && ~isempty(fg))
    warning ('ArgumentCheck:Warning', 'User-supplied gradient and hessResult are not supported for rk4! Finite differences approximation will be used instead.');
end

if ( ~isempty(fg) && ~strcmp(method, 'rk4') )
    options.optOptions = optimset (options.optOptions, 'GradConstr','on');
else
    options.optOptions = optimset (options.optOptions, 'GradConstr','off');
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
        fParam = [];
        if ( options.checkJacobian )
            fParam = f;
        end
        [~, ~, ~, cej_] = ddeParamEst_constraints(fParam, fg, [], tGrid, x, n, np, delays, xHistory, method, options);
    end

    function lsfh = lsf_lh (x, lm_k)
        if ( strcmp(options.hessian_method, 'newton') )
            fParam = [];
            if ( options.checkHessian )
                fParam = f;
            end
            [~, ~, ~, ~, ~, ceh] = ddeParamEst_constraints(fParam, [], fh, tGrid, x, n, np, delays, xHistory, method, options, lm_k(length(x)+1:end));
            % todo serso: check SIGN
            lsfh = H_by_2 + ceh;
        elseif ( strcmp(options.hessian_method, 'gauss-newton') )
            lsfh = H_by_2;
        else
            throw (MException ('ArgumentCheck:IllegalArgument', 'Unsupported method!'));
        end
    end

    function ce_ = ce (x) 
        [~, ce_] = ddeParamEst_constraints(f, [], [], tGrid, x, n, np, delays, xHistory, method, options);
    end

    function [ci, ce, cij, cej] = constraints0 (x) 
        [ci, ce, cij, cej] = ddeParamEst_constraints(f, fg, [], tGrid, x, n, np, delays, xHistory, method, options);
        cej = cej';
    end

if ( ~isempty(fh) && ~strcmp(method, 'rk4') )
    options.optOptions = optimset (options.optOptions, 'Hessian','user-supplied', 'HessFcn', @lsf_lh);
else
    options.optOptions = optimset (options.optOptions, 'Hessian','off');
end



%% solving the task

% sqpTimerId = tic;

if ( options.sqp )
    
    options.sqpOptions.dxmin = options.xTol;
    options.sqpOptions.iterativeTol = options.xTol;
    options.sqpOptions.pSol = options.pSol;
    
    [solution, lmResult, iterations, funCount, ~, sqpInfo] = ddeParamEst_sqp(    ...
        n + np, n - 1, ...
        @(x_k)lsf(x_k), ...
        @(x_k, lm_k)lsf_lh(x_k, lm_k), ...
        @(x_k)ce(x_k), ...
        @(x_k)cej(x_k), ...
        x0, ...
        options.sqpOptions, ...
        debug, ...
        lb, ub);
    
    sosResult = lsf(solution);
    exitflagResult = sqpInfo.flag;
    infoResult = struct('iterations', iterations, 'funcCount', funCount, 'sqpInfo', sqpInfo);
    
else
    
    options.optOptions.TolX = options.xTol / 10;
    
    [solution, sosResult, exitflagResult, infoResult, lmResult] = ...
        fmincon( ...
        @(x_k)lsf(x_k), ...
        x0, [], [], [], [], ...
        lb, ub, ...
        @(x_k)constraints0(x_k), ...
        options.optOptions);
    
end

% sqpTimeResult = toc(sqpTimerId);

% pResult = solution( n + 1 : n + np );
% sqpThetaResult = sqpSolution( n + 1 : n + np );
% 
% display('ThetaResult - SqpThetaResult');
% display(pResult - sqpThetaResult);
% 
% display(sprintf('Time for fmincon: %0.3f s', fminconTimeResult));
% display(sprintf('Time for sqp: %0.3f s', sqpTimeResult));

%% preparing result
pResult = solution( n + 1 : n + np );
% sqpThetaResult = sqpSolution( n + 1 : n + np );

% display('ThetaResult - SqpThetaResult');
% display(pResult - sqpThetaResult);

if ( showResult )
    figure;
    grid on;
    hold on;
    title(sprintf('Method: %s\nApproximation grid: %i\nTime: %0.3f s', method, n, toc(timerId)));
    plot (t, x, '-b');
    plot (utils.interpolate(t, n, 'spline'), solution(1:n), 'or');
end

xResult = zeros(length(t), 1);
j = 1;
for i = 1 : 1 : n
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
    %display(xResult);
    display(sosResult);
end

if ( ~isempty(options.pSol) && isfield(infoResult, 'sqpInfo') && isfield(infoResult.sqpInfo, 'pSolDiffs') )
    oldFormat = get(0,'format');
    format('long');
    display(infoResult.sqpInfo.pSolDiffs');
    format(oldFormat);
end

if (debug || showResult)
    display(sprintf('Method: %s', method));
    display(sprintf('Approximation grid: %i', n));
    
    display(infoResult);
    display(pResult);
end

end