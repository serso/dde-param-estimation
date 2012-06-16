function [x, p, info] = ddeParamEst_TEST ( ...
    nx, ...
    tMin, tMax, ...
    xSolH, pSol, ...
    f, fg, fh, ...
    delays, xHistory, maxDelay, ...
    options, ...
    xSigmaError, tSigmaError, ...
    pLb, pUb, p0)
%DDEPARAMEST_TEST makes simple test of estimating DDE parameters

%% INIT
% setting default values for not obligatory arguments
obligatoryArgs = 14;

np = length(pSol);

if ( nargin <= obligatoryArgs || isempty(pLb) )
    pLb = -Inf * ones (np, 1);
end

if (nargin <= obligatoryArgs + 1 || isempty(pUb))
    pUb = Inf * ones(np, 1);
end

if (nargin <= obligatoryArgs + 2 || isempty(p0))
    p0 = [];
end

%%
% create grid (t, x) and (t, xErr)

[t, xSol, xSolErr, deltaT, ~] = ...
    ddeParamEst_createInitialGrid ( ...
    xSolH, ...
    nx, tMin, tMax, ...
    xSigmaError, tSigmaError);

%% SOLVE
[x, p, info] = ddeParamEst(t, xSolErr, f, fg, fh, delays, xHistory, maxDelay, options, np, pLb, pUb, p0, deltaT);

%% PRINT RESULT
if ( options.showResult )
    display('Known p:');
    display(pSol);
    display('p* (estimated parameters):');
    display(p);
    fprintf('norm(p - p*) = %f\n', norm(pSol - p, inf));
end
% saveas(h, strcat('output/', taskName, '_result'), 'png');

end