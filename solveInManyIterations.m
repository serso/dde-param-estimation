function [x, xResult, thetaResult, sumOfSquares] = solveInManyIterations ( ...
    taskName, N, ...
    tMin, tMax, ...
    exactX, ...
    f, gradF, hessF, ...
    delays, delayF, maxDelay, ...
    method, options, theta, p, ...
    xSigmaError, tSigmaError, ...
    debug, showResult, showIntermidiateResult, thetaLb, thetaUb, x0)

% setting default values for not obligatory arguments
obligatoryArgs = 20;

if ( nargin <= obligatoryArgs || isempty(thetaLb) )
    thetaLb = -Inf * ones (p, 1);
end

if (nargin <= obligatoryArgs + 1 || isempty(thetaUb))
    thetaUb = Inf * ones(p, 1);
end

if (nargin <= obligatoryArgs + 2 || isempty(x0))
    x0 = 0;
end

% create grid (t, x)

[t, xExactSolution, xWithErrors, ~] = ...
    createInitialGrid ( ...
    exactX, ...
    N, tMin, tMax, ...
    xSigmaError, tSigmaError);

[x, xResult, thetaResult, sumOfSquares] = solveInManyIterations1(taskName, t, xWithErrors, f, gradF, hessF, delays, delayF, maxDelay, method, options, p, debug, showResult, showIntermidiateResult, thetaLb, thetaUb, x0);

i = 1;
while ( t (i) - tMin <= maxDelay )
    i = i + 1;
end

h = plot (t(i+1:length(t)), xExactSolution(i+1:length(xExactSolution)), '-r');
saveas(h, strcat('output/', taskName, '_result'), 'png'); 

end