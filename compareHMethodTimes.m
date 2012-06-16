function compareHMethodTimes(tMin, tMax, xSol, f, fg, fh, delays, delayF, maxDelay, options, pSol, xSigmaError, tSigmaError, pLb, pUb, p0, hMethods, ns )

%%

options.showResult = false;
options.plotResult = false;

if ( nargin < 17 || isempty(hMethods) )
    hMethods = options.hessian_method;
end

if ( nargin < 18 || isempty(ns) )
   ns = [100 250 500 750 1000 1500 2000 5000 ];
end

% hMethods = {'fmincon', 'symrcm'};
colors = {'m', 'r', 'g', 'b', 'k', 'c', 'y', 'm'}';
close all;

[p, ~] = size(pSol);
        
figure('Position', [1, 1, 1024, 600]);
grid on;
hold on;

times = Inf * ones(length(ns), length(hMethods));
pErrors = Inf * ones(length(ns), length(hMethods));


for ni = 1 : length(ns)
    
    N = ns(ni);
    
    fprintf('\n###Start step for N');
    fprintf('\n###N = %i', N);
    
    [t, ~, xWithErrors, deltaT, ~] = ...
        ddeParamEst_createInitialGrid ( ...
        xSol, ...
        N, tMin, tMax, ...
        xSigmaError, tSigmaError);


    for methodi = 1 : length(hMethods)
        hMethod = hMethods(methodi);
        
        fprintf('\n######Start step for hMethod');
        display(hMethod);
        
        options.hessian_method = hMethod;
        
        timerId = tic;
        [~, pRes] = ddeParamEst(t, xWithErrors, f, fg, fh, delays, delayF, maxDelay, options, p, pLb, pUb, p0, deltaT);
        times(ni, methodi) = toc(timerId);
        
        pErrors(ni, methodi) = norm(pSol - pRes, inf);
        
        fprintf('\n######End step for hMethod');
    end
    
    fprintf('\n###End step for N');
end

figure;
hold on;
grid on;

times = times';
for methodi = 1 : length(hMethods)
    [colorsLength, ~] = size(colors);
    plot (ns, times(methodi,:), colors{ mod(methodi, colorsLength) + 1});
end

legend(hMethods);

display(pErrors);


%%
end