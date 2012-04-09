function compareTimes(tMin, tMax, x_sol, f, fg, fh, delays, delayF, maxDelay, options, theta, xSigmaError, tSigmaError, thetaLb, thetaUb, theta0, methods, ns )

%%

options.showResult = false;

% ns = [100 250 500 750 1000 1250 1500 1750 2000 2500 3000 5000 ];

if ( nargin < 17 || isempty(methods) )
    methods = {'fmincon', 'default', 'symrcm', 'amd', 'colamd', 'colperm', 'dmperm', 'symamd'};
end

if ( nargin < 18 || isempty(ns) )
   ns = [100 250 500 750 1000 1500 2000 5000 ];
end

% methods = {'fmincon', 'symrcm'};
colors = {'m', 'r', 'g', 'b', 'k', 'c', 'y', 'm'};
close all;
        
figure('Position', [1, 1, 1024, 600]);
grid on;
hold on;
        
for i = 1 : length(methods)
    method = methods(i);
    times = [];
    
    if ( strcmp(method, 'fmincon') ) 
        options.sqp = false;
    else
        options.sqp = true;
        options.sqpOptions.stepMethod = method;
    end
    
    for j = 1 : length(ns)
        if ( strcmp(method, 'fmincon') ) 
           if ( ns(j) >= 1000 ) 
%                times(j) = times(j-1);
%                continue;
                    break;
           end
        end
        [~, ~, ~, ~, times(j)] = ddeParamEst_TEST(...
            ns(j),  tMin, tMax, ...
            x_sol, ...
            f, ...
            fg, ...
            fh, ...
            delays, delayF, maxDelay, ...
            options, theta, length(theta), xSigmaError, tSigmaError, thetaLb, thetaUb, theta0);
    end
    
    plot (ns(1:length(times)), times, colors{i});
end

legend(methods);


%%
end