function [x, xResult, thetaResult, sumOfSquares, time] = ddeParamEst( ...
    t, x, ...
    f, fg, fh, ...
    delays, delayF, maxDelay, ...
    options, p, ...
    thetaLb, thetaUb,theta0, deltaT)

addPath('./sqp/src');
addPath('./reduce');

timerId = tic;

% setting default options parameters 
if ( ~isfield(options, 'debug') )     
    options.debug = false;
end

if ( ~isfield(options, 'showResult') )     
    options.showResult = false;
end

if ( ~isfield(options, 'showIntermidiateResult') )     
    options.showIntermidiateResult = false;
end

if( ~isfield(options, 'method') )
    options.method = 'backward_euler';
end


if ( ~isfield(options, 'sqp') )     
    options.sqp = false;
end    

if ( ~isfield(options, 'xTol') )     
    options.xTol = 10^-2;
end    

if ( ~isfield(options, 'thetaTol') )     
    options.thetaTol = 10^-2;
end    


if ( ~isfield(options, 'maxNumberOfIterations') )     
    options.maxNumberOfIterations = 10;
end    

if ( options.sqp )
    if ( ~isfield(options, 'sqpOptions') )
        sqpOptions.algo_method        = 'Newton';
        sqpOptions.algo_globalization = 'line-search';
        sqpOptions.stepMethod = 'default';
        options.sqpOptions = sqpOptions;
    end
    
    if ( ~isfield(options.sqpOptions, 'algo_method'))
        options.sqpOptions.algo_method        = 'Newton';
    end
    
    if ( ~isfield(options.sqpOptions, 'algo_globalization'))
        options.sqpOptions.algo_globalization        = 'line-search';
    end
    
    if ( ~isfield(options.sqpOptions, 'stepMethod') ) 
        options.stepMethod = 'default';
    end
end

% setting default values for not obligatory arguments
obligatoryArgs = 10;

if ( nargin <= obligatoryArgs || isempty(thetaLb) )
    thetaLb = -Inf * ones (p, 1);
end

if (nargin <= obligatoryArgs + 1 || isempty(thetaUb))
    thetaUb = Inf * ones(p, 1);
end

if (nargin <= obligatoryArgs + 2 || isempty(theta0))
    theta0 = zeros(p, 1);
end

if (nargin <= obligatoryArgs + 3 || isempty(deltaT))
    
    deltaT = Inf;
   
    for i = 1: 1: size(t, 1)
        if ( i > 1 )
            delta = t(i, 1) - t(i - 1, 1);
            if ( delta < deltaT )
                deltaT = delta;
            end           
        end
    end
    
end

if ( options.debug || options.showResult )
    display(strcat('Solving task: ', options.taskName));
    display(strcat('ODE parameter estimation: ', options.method));
end

% create delay function (is not exists)

[t, x, minTDist, delayF] = extractDelayFunctionFromGrid(t, x, delayF, maxDelay);

% +1 has to guarantees that there would not be elements on the bounds
approximationN = ceil( (max(t) - min(t)) / minTDist ) + 1;

if ( options.debug )
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
            
    if ( options.showIntermidiateResult )
        display('Start iteration of algorithm!');
        display(sprintf('Iteration number: %i', i));
        display(sprintf('Number of elements in grid: %i', approximationN));
        display(sprintf('Number of known elements in grid: %i', length(t)));
    end
    
    NGrid(i) = approximationN;
    
    [xResult, thetaResult, sumOfSquares, ~, output, ~, ~, ~, timeResult] = ...
        ddeParamEstStep ( ...
        f, ...
        fg, ...
        fh, ...
        t, x, approximationN, ...
        p, delays, delayF, options.method, options.debug, ...
        options.showIntermidiateResult, options, x0, ...
        deltaT, ...
        thetaLb, thetaUb);
%     
%     absoluteThetaErrors(i) = norm ( theta - thetaResult, inf );
%     absoluteXErrors(i) = norm ( x - xResult, inf );
    iterations(i) = output.iterations;
    funCounts(i) = output.funcCount;
    times(i) = timeResult;
    
    if ( options.showIntermidiateResult )
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
    
    if ( i >= options.maxNumberOfIterations || ...
            (xDiffs(i) < options.xTol && thetaDiffs(i) < options.thetaTol ))
        % time to stop
        break;
    end
    
end

if ( options.showResult && ~options.showIntermidiateResult )
    if (~isempty(xResult))
        figure('Position', [1, 1, 1024, 600]);
        grid on;
        hold on;
        %title(sprintf('Task: %s\nMethod: %s\nApproximation grid: %i\nTime: %0.3f s', options.taskName, options.method, approximationN, toc(timerId)));
        %title(sprintf('Task: %s', options.taskName));
        %plot (t, x, 'xr');
        h = plot (interpolate(t, length(xResult), 'spline'), xResult(1:length(xResult)), '-b');
        xlabel('t, years');
        ylabel('x, billions');
%         saveas(h, strcat('output/', options.taskName, '_result'), 'png'); 
    end
    
    display(output);
end

time = toc(timerId);

if ( options.debug || options.showResult )
    display(strcat('Results for: ', options.taskName));
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
    display(sprintf('Time: %0.5f s', time));
end
end