function compareTimes(tMin, tMax, xSol, f, fg, fh, delays, delayF, maxDelay, options, pSol, xSigmaError, tSigmaError, pLb, pUb, p0, methods, ns )

%%

options.showResult = false;
options.plotResult = false;

% ns = [100 250 500 750 1000 1250 1500 1750 2000 2500 3000 5000 ];

if ( nargin < 17 || isempty(methods) )
    methods = {'fmincon', 'default', 'symrcm', 'amd', 'colamd', 'colperm', 'dmperm', 'symamd'};
end

if ( nargin < 18 || isempty(ns) )
   ns = [100 250 500 750 1000 1500 2000 5000 ];
end

% methods = {'fmincon', 'symrcm'};
colors = {'m', 'r', 'g', 'b', 'k', 'c', 'y', 'm'}';
close all;

[p, ~] = size(pSol);
        
figure('Position', [1, 1, 1024, 600]);
grid on;
hold on;

times = Inf * ones(length(ns), length(methods));
pErrors = Inf * ones(length(ns), length(methods));


iterativeMethods = {'bicgstab','bicgstabl','cgs','gmres', 'lsqr', 'minres', 'qmr', 'symmlq', 'bicg'}';

for ni = 1 : length(ns)
    
    N = ns(ni);
    
    fprintf('\n###Start step for N');
    fprintf('\n###N = %i', N);
    
    [t, ~, xWithErrors, deltaT, ~] = ...
        ddeParamEst_createInitialGrid ( ...
        xSol, ...
        N, tMin, tMax, ...
        xSigmaError, tSigmaError);


    for methodi = 1 : length(methods)
        method = methods(methodi);
        
        fprintf('\n######Start step for method');
        % fprintf('\n######Method = %s', method);
        display(method);
        
        options.sqpOptions.stepMethodIterative = false;
        if ( strcmp(method, 'fmincon') )
            options.sqp = false;
        elseif (ismember(method, iterativeMethods))
            options.sqpOptions.stepMethodIterative = true;
            options.sqp = true;
            options.sqpOptions.stepMethod = method;
        else
            options.sqp = true;
            options.sqpOptions.stepMethod = method;
        end
        
%         if ( strcmp(method, 'fmincon') )
%             if ( ns(ni) >= 5000 )
%                 %                times(j) = times(j-1);
%                 %                continue;
%                 continue;
%             end
%         end
        
        timerId = tic;
        [~, pRes] = ddeParamEst(t, xWithErrors, f, fg, fh, delays, delayF, maxDelay, options, p, pLb, pUb, p0, deltaT);
        times(ni, methodi) = toc(timerId);
        
        pErrors(ni, methodi) = norm(pSol - pRes, inf);
        
        fprintf('\n######End step for method');
    end
    
    fprintf('\n###End step for N');
end

figure;
hold on;
grid on;

times = times';
for methodi = 1 : length(methods)
    [colorsLength, ~] = size(colors);
    plot (ns, times(methodi,:), colors{ mod(methodi, colorsLength) + 1});
end

legend(methods);

display(pErrors);


%%
end