function [x, xResult, thetaResult, sumOfSquares, time] = ddeParamEst_TEST ( ...
    N, ...
    tMin, tMax, ...
    x_sol, ...
    f, fg, fh, ...
    delays, delayF, maxDelay, ...
    options, theta, p, ...
    xSigmaError, tSigmaError, ...
    thetaLb, thetaUb, theta0)

% setting default values for not obligatory arguments
obligatoryArgs = 15;

if ( nargin <= obligatoryArgs || isempty(thetaLb) )
    thetaLb = -Inf * ones (p, 1);
end

if (nargin <= obligatoryArgs + 1 || isempty(thetaUb))
    thetaUb = Inf * ones(p, 1);
end

if (nargin <= obligatoryArgs + 2 || isempty(theta0))
    theta0 = [];
end

% create grid (t, x)

[t, xExactSolution, xWithErrors, deltaT, ~] = ...
    createInitialGrid ( ...
    x_sol, ...
    N, tMin, tMax, ...
    xSigmaError, tSigmaError);

[x, xResult, thetaResult, sumOfSquares, time] = ddeParamEst(t, xWithErrors, f, fg, fh, delays, delayF, maxDelay, options, p, thetaLb, thetaUb, theta0, deltaT);

if (  options.showResult )
    display('Theta*');
    display(theta);
    display('max(theta - thetaResult)');
    display(norm(theta - thetaResult, inf));
end

i = 1;
while ( t (i) - tMin <= maxDelay )
    i = i + 1;
end

if (  options.showResult )
    plot (t(i+1:length(t)), xExactSolution(i+1:length(xExactSolution)), '-r');
    plot (t(i+1:length(t)), xWithErrors(i+1:length(xWithErrors)), '.r');
end
% saveas(h, strcat('output/', taskName, '_result'), 'png'); 

end