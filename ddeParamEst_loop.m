function [xResult, pResult, info] = ddeParamEst_loop(  ...
        f, ...
        fg, ...
        fh, ...
        t, x, n, ...
        np, delays, xHistory, options, ...
        deltaT, ...
        pLb, pUb, p0 )
%DDEPARAMEST_LOOP Function iterates over different approximation densities
%of constraints

%% INIT

% iterations(i) - number of iterations on the i-th step algorithm
iterations = zeros(1, 1);
% funCounts(i) - number of function calls in i-th the step algorithm
funCounts = zeros(1, 1);

% times(i) - time spent on the i-th step
times = zeros(1, 1);

sumsOfSquares = ones(1, 1);


% xDiffs(i) - difference between x*(i-1) and x*(i): norm(x*(i) - x*(i-1), inf)
xDiffs = ones(1, 1);
xDiffs(1) = Inf;

% pDiffs(i) - difference between np*(i-1) and np*(i): norm(np*(i) - np*(i-1),
% inf)
pDiffs = ones(1, 1);
pDiffs(1) = Inf;

% pSolDiffs(i) - difference between p* and np*(i): norm(np*(i) - p*, inf)
pSolDiffs = ones(1, 1);
pDiffs(1) = Inf;

% ns(i) - n on the i-th step
ns = zeros(1, 1);

% init x0 - interpolate given x values
x0 = zeros(n + np, 1);
x0(1 : n ) = utils.interpolate(x, n, 'spline');
x0(n + 1 : n + np) = p0;

prevXResult = [];
prevPResult = [];

pResult = [];
xResult = [];

%% LOOP

sumOfSquares = [];

for i = 1:2147483647
        
    if ( options.showIntermidiateResult )
        fprintf('######### Start step %i\n', i);
        fprintf('Number of elements in approximation grid: %i\n', n);
        fprintf('Number of known elements in grid: %i\n', length(t));
    end
    
    timerId = tic;

    %% DO STEP
    [xResult, pResult, sumOfSquares, ~, output, ~] = ...
        ddeParamEst_step ( ...
        f, ...
        fg, ...
        fh, ...
        t, x, n, ...
        np, delays, xHistory, options.method, options.debug, ...
        options.showIntermidiateResult, options, x0, ...
        deltaT, ...
        pLb, pUb);
    
    %% SAVE STEP STATISTICS

    times(i) = toc(timerId);
    ns(i) = n;
    iterations(i) = output.iterations;
    funCounts(i) = output.funcCount;
    sumsOfSquares(i) = sumOfSquares;
    
    relativeXDiff = Inf;
    if ( ~isempty(prevXResult) )
        xDiffs(i) = norm(prevXResult - xResult, inf);
        relativeXDiff = xDiffs(i) / norm(prevXResult, inf) ;
    end
    
    relativePDiff = Inf;
    if ( ~isempty(prevPResult) )
        pDiffs(i) = norm(prevPResult - pResult, inf);
        relativePDiff = pDiffs(i) / norm(prevPResult, inf);
    end
    
    if ( ~isempty(options.pSol) )
        pSolDiffs(i) = norm(pResult - options.pSol, inf);
    end
    
    %% SIDPLAY INTERMIDIATE DATA
    if ( options.showIntermidiateResult )
        fprintf('\n#########End step %i\n', i');
        
        if ( isfield(output, 'sqpInfo') )
            oldFormat = get(0,'format');
            format('long');

            %fprintf('SQP iteration times (full):\n')
            %display(output.sqpInfo.times);
            fprintf('SQP iteration times (only iteration algorithm):\n')
            display(output.sqpInfo.stepAlgoTimes');
            fprintf('SQP iteration times (only iteration algorithm) MEAN:\n')
            display(mean(output.sqpInfo.stepAlgoTimes));

            format(oldFormat);
        end
    
        display(output);
        display(xResult);
        display(pResult);
        display(xResult);
        display(prevXResult - xResult);
    end
    
    %% CHECK EXIT CONDITIONS
    
    % maximum number of iterations rached
    maxIterationsReached = i >= options.maxNumberOfIterations;
    
    % maximum approximation grid density reached
    maxApproximationGridNReached = n >= options.maxApproximationN;
    
    % tolerance condition
    tolCond = xDiffs(i) < options.xTol && pDiffs(i) < options.pTol;
    
    if ( maxIterationsReached || ...
            tolCond || ...
            maxApproximationGridNReached )
        
        fprintf('\nddeParamEst: stop. Reason: ');
        if ( maxIterationsReached )
            fprintf('Max iterations reached');
        elseif (tolCond)
            fprintf('Tolerance condition');
        elseif (maxApproximationGridNReached)
            fprintf('Max approximation grid density reached');
        end
        
        fprintf('\n            iteration                   = %i     max = %i', i, options.maxNumberOfIterations);
        fprintf('\n            approximation N             = %i     max = %i', n, options.maxApproximationN);
        fprintf('\n            |x_(i-1)-x_(i)|             = %11.5e tol = %11.5e', xDiffs(i), options.xTol);
        fprintf('\n            |p_(i-1)-p_(i)|             = %11.5e tol = %11.5e', pDiffs(i), options.pTol);
        fprintf('\n');
        
        % time to stop
        break;
    end
    
    %% PREPARE NEW STEP
    n = (2 ^ i) * length(t) - 1;
    % n = (2 ^ i) * length(x) - 1;
    
    % interpolate current result to x0 (to initial vector of algorithm)
    x0 = utils.interpolate(xResult, n, 'spline')';
    
    % setting additional initial parameter for p
    x0(n + 1 : n + np) = pResult;
    
    prevPResult = pResult;
    prevXResult = xResult;
    
end

info.sumsOfSquares = sumsOfSquares;
info.sumOfSquares = sumOfSquares;
info.xDiffs = xDiffs;
info.pDiffs = pDiffs;
info.pSolDiffs = pSolDiffs;
    
info.iterations = iterations;
info.funCounts = funCounts;
    
info.ns = ns;
info.times = times;
info.output = output;

end

