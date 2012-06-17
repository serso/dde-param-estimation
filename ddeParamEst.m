function [xResult, pResult, info] = ddeParamEst( ...
    t, x, ...
    f, fg, fh, ...
    delays, xHistory, maxDelay, ...
    options, pn, ...
    pLb, pUb, p0, deltaT)

timerId = tic;

%% ADD PATHS
addpath('./sqp/src');
addpath('./reduce');
addpath('./bd');

%% SET DEFAULT OPTIONS

% setting default options parameters
options = utils.setDefaultOptions(options, ...
    {
    {'debug', false}
    {'showResult', false}
    {'plotResult', false}
    {'plotExtResult', false}
    % applied below
    %{'extTMax', t(end)}
    {'showIntermidiateResult', false}
    {'method', 'backward-euler'}
    {'sqp', true}
    {'xTol', 0.01}
    {'pTol', 0.01}
    {'pSol', []}
    {'checkHessian', false}
    {'checkJacobian', false}
    {'maxApproximationN', 20000}
    {'maxNumberOfIterations', 20}
    {'sqpOptions', []}
    });


if ( options.sqp )
    options.sqpOptions = utils.setDefaultOptions(options.sqpOptions, ...
        {
        {'algo_method', 'Newton'}
        {'algo_globalization', 'line-search'}
        {'stepMethod', 'default'}
        {'stepMethodIterative', false}
        });
end

%% SET DEFAULT FUNCTION ARGUMENTS

% setting default values for not obligatory arguments
obligatoryArgs = 10;

if ( nargin <= obligatoryArgs || isempty(pLb) )
    pLb = -Inf * ones (pn, 1);
end

if (nargin <= obligatoryArgs + 1 || isempty(pUb))
    pUb = Inf * ones(pn, 1);
end

if (nargin <= obligatoryArgs + 2 || isempty(p0))
    p0 = zeros(pn, 1);
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
%% SHOW ALGORITHM INFORMATION

fprintf('\n######Starting ddeParamEst with next initial parameters: ');
display(options);
if ( options.sqp && isfield(options, 'sqpOptions') )
    display(options.sqpOptions);
end
if ( ~options.sqp && isfield(options, 'optOptions') )
    display(options.optOptions);
end

%% SOLVE

% create delay function (if not exists)

origT = t;
origX = x;
[t, x, minTDist, xHistory] = ddeParamEst_extractXHistoryFromGrid(t, x, xHistory, maxDelay);

% t may be recalculated => apply extTMin and extTMax only here
options = utils.setDefaultOptions(options, ...
    {
    {'extTMax', t(end)}
    });

% +1 has to guarantees that there would not be elements on the bounds
approximationN = ceil( (max(t) - min(t)) / minTDist ) + 1;

% if ( options.debug )
%     display(x);
%     display(t);
% end

[xResult, pResult, info] = ddeParamEst_loop( ...
    f, ...
    fg, ...
    fh, ...
    t, x, approximationN, ...
    pn, delays, xHistory, options, ...
    deltaT, ...
    pLb, pUb, p0);

time = toc(timerId);

%% SHOW RESULT
if ( options.plotResult )
    
    if (~isempty(xResult))
        
        if ( isfield(options, 'extFigureHandle') && ~isempty(options.extFigureHandle) )
            axes(options.extFigureHandle);
            cla(gca);
            grid on;
            hold on;
        else
            % plot result
            figure('Position', [1, 1, 1024, 600]);
            grid on;
            hold on;
            
            h = gca;
            set(h, 'FontSize', 18);
            
            xlabel('t', 'FontSize', 18);
            ylabel('x', 'FontSize', 18);
        
        end;
        
        %title(sprintf('Task: %s\nMethod: %s\nApproximation grid: %i\nTime: %0.3f s', options.taskName, options.method, approximationN, toc(timerId)));
        %title(sprintf('Task: %s', options.taskName));
        %plot (t, x, 'xr');
        
        tResult = utils.interpolate(t, length(xResult), 'spline');
        resultDataPlot = ddeParamEst_plot(tResult, xResult, '-b', options);
        inputDataPlot = ddeParamEst_plot(origT, origX, 'xr', options);
        
        [~, ~, ~, plotTextStrings] = legend('DDE parameter estimation result', 'Input', 'Location', 'Best');
        % title(options.taskName, 'FontSize', 18);
        %         saveas(h, strcat('output/', options.taskName, '_result'), 'png');
        
        
        % plot extrapolation to the right
        if ( options.plotExtResult )
            
            tSpan = [t(1), options.extTMax];
            
            if ( isempty(delays) || ( length(delays) == 1 && delays(1) == 0 ) )
                odeOptOptions = odeset('AbsTol', min([options.xTol, options.pTol]));
                
                odeF = @(t, x, p) f(x', t, p);
                
                xdeResult = ode45(odeF, tSpan, x(1), odeOptOptions, pResult);
            else
                
                % prepare DDE input data
                ddeOptOptions = ddeset('AbsTol', min([options.xTol, options.pTol]));
                
                % convert ddeParamEst input to dde23 input
                ddeF = @(t, x, delays, p) f([x delays]', t, p);
                
                if (length(delays) == 1)
                    ddeDelays = delays(1);
                else
                    ddeDelays = delays(2:end);
                end
                
                ddeXHistory = @(t, p) xHistory(t);
                
                % do extrapolation
                xdeResult = dde23(ddeF, ddeDelays, ddeXHistory, tSpan, ddeOptOptions, pResult);
            end
            
            extT = linspace(tSpan(1), tSpan(2), 1000);
            extX = deval(xdeResult, extT);
            
            
            extDataPlot = ddeParamEst_plot(extT, extX, '-g', options);
            legend([resultDataPlot; inputDataPlot; extDataPlot] , {plotTextStrings{:}, 'Extrapolated result'});
        end
    end
end

if ( options.debug || options.showResult )
    
    display('######################################');
    display(strcat('Results for: ', options.taskName));
    fprintf('Approximation grid: %i\n', approximationN);
    
    fprintf('Last output of step algorithm\n');
    disp(info.output);
    
    %     display(absoluteXErrors);
    %     display(absoluteThetaErrors);
    
    oldFormat = get(0,'format');
    format('long');
    
    if (~isempty(pResult))
        disp('P*: ');
        disp(pResult);
    end
    
    if (~isempty(info.sumOfSquares))
        fprintf('Sum of squares: %d\n', info.sumOfSquares);
    end
    
    fprintf('Sum of squares on the i-th step\n');
    disp(info.sumsOfSquares);
    
    fprintf('X differences: |x_(i-1)-x_(i)| where i is i-th step of algorithm\n');
    disp(info.xDiffs);
    
    fprintf('P differences: |p_(i-1)-p_(i)| where i is i-th step of algorithm\n');
    disp(info.pDiffs);
    
    if ( ~isempty(options.pSol) )
        fprintf('P solution differences: |p*-p_(i)| where i is i-th step of algorithm\n');
        disp(info.pSolDiffs);
    end
    
    format(oldFormat);
    
    fprintf('Iterations of step algorithm\n');
    disp(info.iterations);
    
    fprintf('Number of times function is called on each step of algorithm\n');
    disp(info.funCounts);
    
    fprintf('Approximation grid density: number of elements in approximation grid per step\n');
    disp(info.ns);
    
    fprintf('Step times\n');
    disp(info.times);
    
    fprintf('Total time: %0.5f s\n', time);
    display('######################################');
end
end