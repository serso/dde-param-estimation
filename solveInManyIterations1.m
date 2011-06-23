function [x, xResult, thetaResult, sumOfSquares] = solveInManyIterations1( ...
    taskName, ...
    t, x, ...
    f, gradF, hessF, ...
    delays, delayF, maxDelay, ...
    method, options, p, ...
    debug, showResult, showIntermidiateResult, thetaLb, thetaUb, x0, theta0)

timerId = tic;

% setting default values for not obligatory arguments
obligatoryArgs = 15;

if ( nargin <= obligatoryArgs || isempty(thetaLb) )
    thetaLb = -Inf * ones (p, 1);
end

if (nargin <= obligatoryArgs + 1 || isempty(thetaUb))
    thetaUb = Inf * ones(p, 1);
end

if (nargin <= obligatoryArgs + 2 || isempty(x0))
    x0 = 0;
end

if (nargin <= obligatoryArgs + 3 || isempty(theta0))
    theta0 = 0;
end

if ( debug || showResult )
    display(strcat('Solving task: ', taskName));
    display(strcat('ODE parameter estimation: ', method));
end

% best candidates for function arguments
xTol = 10^-2;
thetaTol = 10^-2;
maxNumberOfIterations = 1;

% create delay function (is not exists)

[t, x, minTDist, delayF] = extractDelayFunctionFromGrid(t, x, delayF, maxDelay);

% +1 has to guarantees that there would not be elements on the bounds
approximationN = ceil( (max(t) - min(t)) / minTDist ) + 1;

if ( debug )
    display(x);
    display(t);
end

% absoluteThetaErrors = ones(1, 1);
% absoluteXErrors = ones(1, 1);

iterations = zeros(1, 1);
funCounts = zeros(1, 1);

times = zeros(1, 1);

xDiffs = ones(1, 1);
xDiffs(1) = Inf;

thetaDiffs = ones(1, 1);
thetaDiffs(1) = Inf;

NGrid = zeros(1, 1);

x0 = zeros(approximationN + p, 1);
x0(1 : approximationN ) = interpolate(x, approximationN, 'spline');
x0(approximationN + 1 : approximationN + p) = theta0;

prevXResult = [];
prevThetaResult = [];

thetaResult = [];
xResult = [];

sumOfSquares = [];

for i = 1:2147483647
    
    if ( showIntermidiateResult )
        display('Start iteration of algorithm!');
        display(sprintf('Iteration number: %i', i));
        display(sprintf('Number of elements in grid: %i', approximationN));
        display(sprintf('Number of known elements in grid: %i', length(t)));
    end
    
    NGrid(i) = approximationN;
    
    [xResult, thetaResult, sumOfSquares, ~, output, ~, ~, ~, timeResult] = ...
        solve ( ...
        f, ...
        gradF, ...
        hessF, ...
        t, x, approximationN, ...
        p, delays, delayF, method, debug, ...
        showIntermidiateResult, options, x0, ...
        thetaLb, thetaUb);
%     
%     absoluteThetaErrors(i) = norm ( theta - thetaResult, inf );
%     absoluteXErrors(i) = norm ( x - xResult, inf );
    iterations(i) = output.iterations;
    funCounts(i) = output.funcCount;
    times(i) = timeResult;
    
    if ( showIntermidiateResult )
        display('Iteration ended!');
        display(sprintf('Iteration number: %i', i));
        display(output);
        display(xResult);
        display(thetaResult);
        display(xResult);
        display(prevXResult - xResult);
    end
    
    if ( ~isempty(prevXResult) )
        xDiffs(i) = norm(prevXResult - xResult, inf);
    end
    
    if ( ~isempty(prevThetaResult) )
        thetaDiffs(i) = norm(prevThetaResult - thetaResult, inf);
    end
    
    prevThetaResult = thetaResult;
    prevXResult = xResult;
    
    %% new dimension of vector x
    approximationN = ( i + 1 ) * length(x) - 1;
    
    % interpolate current result to x0 (to initial vector of algorithm)
    x0 = interpolate(xResult, approximationN, 'spline')';
    
    % setting additional initial parameter for theta
    x0(approximationN + 1 : approximationN + p) = thetaResult;
    
    if ( i >= maxNumberOfIterations || ...
            xDiffs(i) < xTol || ...
            thetaDiffs(i) < thetaTol )
        % time to stop
        break;
    end
    
end

if ( showResult && ~showIntermidiateResult )
    if (~isempty(xResult))
        figure('Position', [1, 1, 1024, 600]);
        grid on;
        hold on;
        %title(sprintf('Task: %s\nMethod: %s\nApproximation grid: %i\nTime: %0.3f s', taskName, method, approximationN, toc(timerId)));
        %title(sprintf('Task: %s', taskName));
        %plot (t, x, 'xr');
        h = plot (interpolate(t, length(xResult), 'spline'), xResult(1:length(xResult)), '-b');
        xlabel('t, years');
        ylabel('x, billions');
%         saveas(h, strcat('output/', taskName, '_result'), 'png'); 
    end
    
    display(output);
end

if ( debug || showResult )
    display(strcat('Results for: ', taskName));
    display(sprintf('Approximation grid: %i', approximationN));
    
%     display(absoluteXErrors);
%     display(absoluteThetaErrors);
    
     if (~isempty(thetaResult))
         display('Theta result');
         display(thetaResult);
     end
    
    if (~isempty(sumOfSquares))
        display(sprintf('Sum of squares: %s', sumOfSquares));
    end
    
    display(xDiffs);
    display(thetaDiffs);
    
    display(iterations);
    display(funCounts);
    
    display(NGrid);
    display(times);
    
%     figure;
%     grid on;
%     hold on;
%     title('Time(N)');
%     plot (NGrid, times, '-b');
%     
    display(sprintf('Time: %0.5f s', toc(timerId)));
end
end