 
function ddeParamEst_TEST_SUITE()

%% intial clear data

clc;
close all;

%% Parameter estimation in ODE

clear;

showResult = true;
clearData = true;

% methods = {'euler' 'backward-euler' 'box' 'rk4'};

% number of known points (i.e. values of function x(t))
N = 200;

xSigmaError = 0.01;
tSigmaError = 0.0;

options.xTol = 0.01;
options.pTol = 0.01;

%optOptions = optimset('Algorithm', 'sqp');
optOptions = optimset('Algorithm', 'interior-point', 'MaxFunEvals', 6000);
% optOptions = optimset(optOptions, 'SubproblemAlgorithm', 'cg');
optOptions = optimset(optOptions, 'DerivativeCheck', 'off');
optOptions = optimset(optOptions, 'FinDiffType', 'central');
%optOptions = optimset(optOptions, 'Algorithm', 'sqp');

options.optOptions = optOptions;

options.sqp = true;
sqpOptions.miter = 30;

options.hessian_method = 'gauss-newton';
%options.hessian_method = 'newton';

options.maxNumberOfIterations = 10;

%sqpOptions.algo_method        = 'quasi-Newton';
sqpOptions.algo_method        = 'Newton';

sqpOptions.algo_globalization = 'line-search';
%sqpOptions.algo_globalization = 'unit step-size';

sqpOptions.stepMethodIterative = false;
sqpOptions.stepMethod = 'ldls';
%sqpOptions.stepMethod = 'bicg';
%sqpOptions.stepMethod = 'block-decomposition';
sqpOptions.ldlsThreshold = 0.001;
sqpOptions.iterativeMaxit = 1000;
sqpOptions.iterativePrecondAlgorithm = 'luinc';
%sqpOptions.iterativePrecondAlgorithm = 'no';
sqpOptions.iterativePrecondAlgorithmThresh = 0.01;
%sqpOptions.iterativePrecondAlgorithmThresh = 0.1;
options.sqpOptions = sqpOptions;

options.checkHessian = false;
options.checkJacobian = false;

%     sqpOptions.tol(1)  = tolopt(1);  % tolerance on the gradient of the Lagrangian
%     sqpOptions.tol(2)  = tolopt(2);  % tolerance on the feasibility
%     sqpOptions.tol(3)  = tolopt(3);  % tolerance on the complementarity

% options.sqpOptions = sqpOptions;

options.method = 'backward-euler';
%options.method = 'box';
%options.method = 'euler';

%options.method = method;
options.debug = false; 
options.showResult = showResult; 
options.plotResult = showResult;
options.plotExtResult = false;

%% intial clear data

clc;
close all;
%% 1

% profile on;

optionsCopy = struct(options);
optionsCopy.taskName = 'task-01';    

if ( clearData )
    clc;    
    close all;
end

tMin = 0;
tMax = 10;

pSol = 2;
optionsCopy.pSol = pSol;

xSolH = @(t) pSol - (pSol - 1) * exp(-t);

f = @(x, ~, p) p - x;
fg = @(~, ~, ~) [ -1, 1 ];
fh = @(x, t, p) [ 0, 0; 0 0];

ddeParamEst_TEST(...
    N,  tMin, tMax, ...
    xSolH, pSol, ...
    f, ...
    fg, ...
    fh, ...
    0, [], [], ...
    optionsCopy, xSigmaError, tSigmaError);

% profile viewer;

%%

%{'fmincon', 'default', 'symrcm', 'amd', 'colamd', 'colperm', 'dmperm',
%'symamd'}

compareTimes(tMin, tMax, ...
    xSolH, ...
    f, ...
    fg, ...
    fh, ...
    0, [], [], ...
    options, pSol,  xSigmaError, tSigmaError, [], [], [], ...
    {'fmincon', 'ldl', 'mldivide', 'bicg', 'cgs'} , ...
    [100 500 1000]);

%%
compareTimes(tMin, tMax, ...
    xSolH, ...
    f, ...
    fg, ...
    fh, ...
    0, [], [], ...
    options, pSol,  xSigmaError, tSigmaError, [], [], [], ...
    {'ldl'} , ...
    [50 100 250 500 750 1000 2500 4000 6000]);

%% 2

options.taskName = 'task_02';

if ( clearData )
    clc;
    close all;
end

tMin = 0;
tMax = 6;

pSol = 2;      

xSolH = @(t) pSol - sin( t ) * exp(-t);

f = @(x, t, p) p - x - exp( - t ) * cos (t);

fg = @(x, t, p) [ -1, 1 ];
fh = @(x, t, p) [ 0, 0; 0 0];

ddeParamEst_TEST(...
    N,  tMin, tMax, ...
    xSolH, pSol, ...
    f, ...
    fg, ...
    fh, ...
    0, [], [], ...
    options, xSigmaError, tSigmaError);

%% 3

options.taskName = 'task_03';

if ( clearData )
    clc;
    close all;
end

tMin = 0.5;
tMax = 10;

pSol = 2;

xSolH = @(t) pSol * cos(t) / t;

f = @(x, t, p) -( p * sin (t) / t + x / t  );

fg = @(x, t, p) [ -1/t, -sin(t)/t ];
fh = @(x, t, p) [ 0, 0; 0, 0];

ddeParamEst_TEST(...
    N,  tMin, tMax, ...
    xSolH, pSol, ...
    f, ...
    fg, ...
    fh, ...
    0, [], [], ...
    options, xSigmaError, tSigmaError);

%% 4

options.taskName = 'task_04';

if ( clearData )
    clc;
    close all;
end

tMin = 0;
tMax = 10;

pSol = [ 4; 2];

xSolH = @(t) pSol(1) / pSol(2) * ( 1 - exp( - pSol(2) * t ));

f = @(x, t, p) p(1) - p(2) * x;

fg = @(x, t, p) [ -p(2), 1, -x ];
fh = @(x, t, p) [
    0,  0, -1;
    0, 0, 0;
    -1, 0, 0];

ddeParamEst_TEST(...
    N,  tMin, tMax, ...
    xSolH, pSol, ...
    f, ...
    fg, ...
    fh, ...
    0, [], [], ...
    options, xSigmaError, tSigmaError);

%% 5

%profile on;

options.taskName = 'task_05';

if ( clearData )
    clc;
    close all;
end

tMin = 0;
tMax = 10;

pSol = [ 2; 1; 1; 1];
p0 = [ 2.1; 0.9; 0.9; 1.1];

pLb = [];%[ 0 0 0 0 ];
pUb = [];%[ 8 5 5 5 ];
%p = 4;
%pUb = Inf * ones(p, 1);
%pLb = - pUb;

xSolH = @(t) pSol(1) - pSol(2) * sin( pSol(3) * t ) * exp (- pSol(4) * t);

%f = @(x, t, pSol) pSol(4) * (pSol(1) - x) - pSol(2) * pSol(3) * cos(pSol(3) * t ) * exp( - pSol(4) * t );

[f, fg, fh] = ddeParamEst_checkUserInput('p_4 * (p_1 - x) - p_2 * p_3 * cos( p_3 * t ) * exp( - p_4 * t )', length(pSol), [], []);

% fg0 = @(x, t, theta_1, theta_2, theta_3, theta_4) [ -theta_4, theta_4, -(theta_3*cos(t*theta_3))/exp(t*theta_4), -(theta_3*cos(t*theta_3))/exp(t*theta_4), theta_1 - x + (t*theta_2*theta_3*cos(t*theta_3))/exp(t*theta_4)]';
% fg = @(x, t, pSol) fg0(x, t, pSol(1), pSol(2), pSol(3), pSol(4)); 
% 
% fh0  = @(x, t, theta_1, theta_2, theta_3, theta_4) ...
%       [  0, 0,                                         0,                                         0,                                                   -1;
%          0, 0,                                         0,                                         0,                                                    1;
%          0, 0,                                         0,                                         0,            (t*theta_3*cos(t*theta_3))/exp(t*theta_4);
%          0, 0,                                         0,                                         0,            (t*theta_3*cos(t*theta_3))/exp(t*theta_4);
%         -1, 1, (t*theta_3*cos(t*theta_3))/exp(t*theta_4), (t*theta_3*cos(t*theta_3))/exp(t*theta_4), -(t^2*theta_2*theta_3*cos(t*theta_3))/exp(t*theta_4)];
% fh = @(x, t, pSol) fh0(x, t, pSol(1), pSol(2), pSol(3), pSol(4)); 
%     
ddeParamEst_TEST(...
    N,  tMin, tMax, ...
    xSolH, pSol, ...
    f, ...
    fg, ...
    fh, ...
    0, [], [], ...
    options, xSigmaError, tSigmaError, pLb, pUb, p0);

%profile viewer;
%%

clc;
compareTimes(tMin, tMax, ...
    xSolH, ...
    f, ...
    fg, ...
    fh, ...
    0, [], [], ...
    options, pSol,  xSigmaError, tSigmaError, [], [], [], ...
    {'ldl'} , ...
    [1000]);

%% 5.5

optionsCopy = struct(options);
optionsCopy.taskName = 'task_05.5';
optionsCopy.xTol = 10^-3;
optionsCopy.plotExtResult = false;

optionsCopy.hessian_method = 'newton';
optionsCopy.method = 'backward-euler';
%optionsCopy.method = 'euler';

if ( clearData )
    clc;
    close all;
end

pSol = 1;
optionsCopy.pSol = pSol;

pUb = [];
pLb = [];
p0 = 0;
%pUb = Inf * ones(p, 1);
%pLb = - pUb;

tMin = 0;
tMax = 1;

odeOptOptions = odeset('AbsTol', min([optionsCopy.xTol, optionsCopy.pTol]));
odef = @(t, x, p) exp(x*p);
odeSol = ode45(odef, [tMin, tMax], 0, odeOptOptions, pSol);
xSolH = @(t) interp1(odeSol.x, odeSol.y, t, 'spline');

[ f, fg, fh, fSym, fgSym fhSym] = ddeParamest_checkUserInput('exp(x*p_1)', 1, 0, []);
pretty(fSym);
pretty(fgSym);
pretty(fhSym);

ddeParamEst_TEST(...
    N,  tMin, tMax, ...
    xSolH, pSol,...
    f, ...
    fg, ...
    fh, ...
    0, xSolH, [], ...
    optionsCopy, 0, tSigmaError, pLb, pUb, p0);

%% 6

% profile on;

options.showResult = true;
options.taskName = 'task_06';

if ( clearData )
    clc;
    close all;
end

pSol = 1;

tau = 2;

tMin = -tau;
tMax = 10;

pUb = [];
pLb = [];

xSolH = @(t) cos( pSol * pi * t / ( 2 * tau));

%f = @(x, t, p) ( - p * pi / (2 * tau)) * x(2);

xHistory = @(t) xSolH(t);

[f, fg, fh] = ddeParamEst_checkUserInput('( - p_1 * pi / (2 * tau)) * x_1', length(pSol), 1, {{'tau', tau}});

% fg = @(x, t, p) [0, - p * pi / (2 * tau), (- pi / (2 * tau)) * x(2)];
% 
% fh = @(x, t, p) [
%     0,                  - pi / (2 * tau); ...
%     - pi / (2 * tau),     0];


delays = [ 0 tau ];

ddeParamEst_TEST(...
    N,  tMin, tMax, ...
    xSolH, pSol, ...
    f, ...
    fg, ...
    fh, ...
    delays, xHistory, max(delays), ...
    options, xSigmaError, tSigmaError, pLb, pUb);

%profile viewer;
%%
compareTimes(tMin, tMax, ...
    xSolH, ...
    f, ...
    fg, ...
    fh, ...
    delays, xHistory, max(delays), ...
    options, pSol,  xSigmaError, tSigmaError, [], [], [], ...
    {'ldl', 'fmincon', 'bicg'} , ...
    [50 100 500 1000 2000 3000 5000 10000]);
%%

clc;
optionsCopy = struct(options);
optionsCopy.plotResult = false;
optionsCopy.xTol = 0.01;
optionsCopy.pTol = 0.01;
compareTimes(tMin, tMax, ...
    xSolH, ...
    f, ...
    fg, ...
    fh, ...
    delays, xHistory, max(delays), ...
    optionsCopy, pSol,  xSigmaError, tSigmaError, [], [], [], ...
    {'ldl', 'bicg'} , ...
    [50 100 250 500 750 1000 2500 4000 6000]);

%%

clc;

compareTimes(tMin, tMax, ...
    xSolH, ...
    f, ...
    fg, ...
    fh, ...
    0, [], [], ...
    options, pSol,  xSigmaError, tSigmaError, [], [], [], ...
    {'ldl'} , ...
    [1000]);

%% 7
options.taskName = 'task_07';

if ( clearData )
    clc;
    close all;
end

pSol = 1;
k = 1;
tau = pi / 2;
b = 1;

tMin = -tau;
tMax = 10;


%pLb = 0;
%pUb = 4;

xSolH = @(t) b * sin ( k * pSol * t ) - b * sin ( k * pSol * tau ) / ( 1 - cos ( k * pSol * tau )) * cos ( k * pSol * t );

xHistory = @(t) xSolH(t);

[f, fg, fh, fSym, fgSym, fhSym] = ddeParamEst_checkUserInput('- p_1 * k * x_1', length(pSol), 1, {{'k', k}});
pretty(fSym);
pretty(fgSym);
pretty(fhSym);

delays = [ 0 tau ];

ddeParamEst_TEST(...
    N,  tMin, tMax, ...
    xSolH, pSol,...
    f, ...
    fg, ...
    fh, ...
    delays, xHistory, max(delays), ...
    options, xSigmaError, tSigmaError);
%%
clc;
close all;

optionsCopy = struct(options);
optionsCopy.maxNumberOfIterations = 1;
optionsCopy.plotResult = false;

optionsCopy.sqpOptions = struct(options.sqpOptions);
optionsCopy.sqpOptions.iterativeMaxit = 1000;
optionsCopy.sqpOptions.stepMethodIterative = true;
optionsCopy.sqpOptions.stepMethod = 'bicg';

[t, ~, xWithErrors, deltaT, ~] = ...
    ddeParamEst_createInitialGrid ( ...
    xSolH, ...
    1000, tMin, tMax, ...
    xSigmaError, tSigmaError);

for iterativePrecondAlgorithmThresh = [0.0001 0.001 0.01 0.02 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2]%[0.0001 0.001 0.01 0.02 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.25 0.3 0.35 0.4 0.5]
    
    optionsCopy.sqpOptions.iterativePrecondAlgorithmThresh = iterativePrecondAlgorithmThresh;
    fprintf('### luinc tolerance = %d', iterativePrecondAlgorithmThresh);
    
    ddeParamEst(t, xWithErrors, f, fg, fh, delays, xHistory, max(delays), optionsCopy, length(pSol), [], [], [], deltaT);    
end

%%
clc;
close all;

N = 1000;
options.maxNumberOfIterations = 1;
options.sqpOptions.iterativeMaxit = 1000;

[t, ~, xWithErrors, deltaT, ~] = ...
    ddeParamEst_createInitialGrid ( ...
    xSolH, ...
    N, tMin, tMax, ...
    xSigmaError, tSigmaError);

for iterativePrecondAlgorithmThresh = [0.0001 0.001 0.01 0.02 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2]%[0.0001 0.001 0.01 0.02 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.25 0.3 0.35 0.4 0.5]
    
    options.sqpOptions.iterativePrecondAlgorithmThresh = iterativePrecondAlgorithmThresh;
    fprintf('### luinc tolerance = %d', iterativePrecondAlgorithmThresh);
    
    ddeParamEst(t, xWithErrors, f, fg, fh, delays, xHistory, max(delays), options, length(pSol), pLb, pUb, [], deltaT);
    close all;
    
end
options.maxNumberOfIterations = 20;
%% 7.5

optionsCopy = struct(options);
optionsCopy.taskName = 'task_07.5';
optionsCopy.xTol = 10^-3;
optionsCopy.plotExtResult = false;

%options.hessian_method = 'gauss-newton';
optionsCopy.hessian_method = 'newton';
optionsCopy.checkJacobian = false;
optionsCopy.checkHessian = true;
optionsCopy.method = 'backward-euler';
%optionsCopy.maxNumberOfIterations = 4;
%options.sqpOptions.ldlsThreshold = 0.000001;

if ( clearData )
    clc;
    close all;
end

pSol = 1;
optionsCopy.pSol = pSol;

pUb = [];
pLb = [];
p0 = 1;
%pUb = Inf * ones(p, 1);
%pLb = - pUb;

delays = [0 0.5 1.9];
ddeDelays = delays(2:end);

tMin = -max(delays);
tMax = 1;

ddef = @(t,y,Z,lambda) sin(y*Z(1)*Z(2)*lambda);
ddeXHistory = @(t,lambda) t;
ddeSol = dde23(ddef,ddeDelays,ddeXHistory,[tMin, tMax],[],pSol);

xSolH = @(t) interp1(ddeSol.x, ddeSol.y, t, 'spline');

[ f, fg, fh, fSym, fgSym fhSym] = ddeParamest_checkUserInput('sin(x*x_1*x_2*p_1)', 1, length(delays), []);
pretty(fSym);
pretty(fgSym);
pretty(fhSym);

ddeParamEst_TEST(...
    N,  tMin, tMax, ...
    xSolH, pSol,...
    f, ...
    fg, ...
    fh, ...
    delays, xSolH, max(delays), ...
    optionsCopy, 0, tSigmaError, pLb, pUb, p0);

%% 8

optionsCopy = struct(options);
optionsCopy.taskName = 'task_08';
optionsCopy.xTol = 10^-1;
%optionsCopy.maxNumberOfIterations = 4;
%options.sqpOptions.ldlsThreshold = 0.000001;

if ( clearData )
    clc;
    close all;
end

pSol = 2;
optionsCopy.pSol = pSol;

pUb = [];
pLb = [];
p0 = 0.5;
%pUb = Inf * ones(p, 1);
%pLb = - pUb;

tau = 1;

tMin = tau;
tMax = 20;

delays = [ 0 tau ];

ddef = @(t,y,Z,lambda) -lambda*Z*(1 + y);
ddeXHistory = @(t,lambda) t;
ddeSol = dde23(ddef,tau,ddeXHistory,[tMin - tau, tMax],[],pSol);

xSolH = @(t) interp1(ddeSol.x, ddeSol.y, t, 'spline');
%xSolH = @(t) interp1q(exampleSolution.x', exampleSolution.y', t);

% f = @(x, t, pSol) -pSol * x(2) * (1 + x(1));
% 
% xHistory = []; %@(t) t;
% 
% fg = @(x, t, pSol) [- pSol * x(2), -pSol * ( 1 + x(1)), -x(2) * (1 + x(1))];
% 
% fh = @(x, t, pSol) [
%     - pSol * 2,   -1 - 2 * x(2); ...
%     -1 - 2 * x(2),           0];
[ f, fg, fh, fSym, fgSym fhSym] = ddeParamest_checkUserInput('-p_1 *x_1 * (1 + x)', 1, 1, []);
pretty(fSym);
pretty(fgSym);
pretty(fhSym);
xHistory = [];

ddeParamEst_TEST(...
    N,  tMin, tMax, ...
    xSolH, pSol,...
    f, ...
    fg, ...
    fh, ...
    delays, xSolH, max(delays), ...
    optionsCopy, xSigmaError, tSigmaError, pLb, pUb, p0);
%%
compareHMethodTimes(tMin, tMax, ...
    xSolH, ...
    f, ...
    fg, ...
    fh, ...
    delays, xHistory, max(delays), ...
    options, pSol,  xSigmaError, tSigmaError, pLb, pUb, p0, ...
    {'newton', 'gauss-newton'} , ...
    [300 450 800]);
%%
clc;
close all;

n = 1000;
optionsCopy = struct(options);
optionsCopy.sqpOptions = struct(options.sqpOptions);
optionsCopy.maxNumberOfIterations = 3;
optionsCopy.xTol = 0.01;
optionsCopy.pTol = 0.01;
optionsCopy.sqpOptions.stepMethodIterative = false;
optionsCopy.sqpOptions.stepMethod = 'ldls';

[t, ~, xWithErrors, deltaT, ~] = ...
    ddeParamEst_createInitialGrid ( ...
    xSolH, ...
    n, tMin, tMax, ...
    xSigmaError, tSigmaError);

times = [];
ps = [];
ldlsThresholds = [0.000000001 0.00000001 0.0000001 0.000001 0.00001 0.0001 0.001 0.01 0.1 0.2 0.3 0.4 0.5];
i = 1;
for ldlsThreshold = ldlsThresholds
    
    timerId = tic;
    optionsCopy.sqpOptions.ldlsThreshold = ldlsThreshold;
    fprintf('### LDL threshold = %d', ldlsThreshold);
    
    [~, ps(i), ~] = ddeParamEst(t, xWithErrors, f, fg, fh, delays, xHistory, max(delays), optionsCopy, length(pSol), pLb, pUb, [], deltaT);
    close all;
    
    times(i) = toc(timerId);
    i = i + 1;
    
end
display(ldlsThresholds');
display(times');
display(ps');

%%
clc;
close all;

N = 1000;
options.maxNumberOfIterations = 3;

[t, ~, xWithErrors, deltaT, ~] = ...
    ddeParamEst_createInitialGrid ( ...
    xSolH, ...
    N, tMin, tMax, ...
    xSigmaError, tSigmaError);

times = [];
thetas = [];
maxits = 3:1:20;
i = 1;
for iterativeMaxit = maxits
    
    timerId = tic;
    options.sqpOptions.iterativeMaxit = iterativeMaxit;
    fprintf('### max iterations in SQP step = %i', iterativeMaxit);
    
    [~, ~, thetas(i), ~, ~] = ddeParamEst(t, xWithErrors, f, fg, fh, delays, xHistory, max(delays), options, length(pSol), pLb, pUb, [], deltaT);
    close all;
    
    times(i) = toc(timerId);
    i = i + 1;
    
end
display(maxits');
display(times');
display(thetas');
options.maxNumberOfIterations = 20;

%%
clc;
close all;

N = 1000;
options.maxNumberOfIterations = 1;
options.sqpOptions.iterativeMaxit = 1000;

[t, ~, xWithErrors, deltaT, ~] = ...
    ddeParamEst_createInitialGrid ( ...
    xSolH, ...
    N, tMin, tMax, ...
    xSigmaError, tSigmaError);

for iterativePrecondAlgorithmThresh = [0.0001 0.001 0.01 0.02 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2]
    
    options.sqpOptions.iterativePrecondAlgorithmThresh = iterativePrecondAlgorithmThresh;
    fprintf('### luinc tolerance = %d', iterativePrecondAlgorithmThresh);
    
    ddeParamEst(t, xWithErrors, f, fg, fh, delays, xHistory, max(delays), options, length(pSol), pLb, pUb, [], deltaT);
    close all;
    
end
options.maxNumberOfIterations = 20;

%%
compareTimes(tMin, tMax, ...
    xSolH, ...
    f, ...
    fg, ...
    fh, ...
    delays, xHistory, max(delays), ...
    options, pSol,  xSigmaError, tSigmaError, [], [], [], ...
    {'ldl', 'fmincon', 'bicg', 'cgs'} , ...
    [50 100 500 1000 5000]);
%%

clc;

compareTimes(tMin, tMax, ...
    xSolH, ...
    f, ...
    fg, ...
    fh, ...
    0, [], [], ...
    options, pSol,  xSigmaError, tSigmaError, [], [], [], ...
    {'ldl'} , ...
    [50 100 250 500 750 1000 2500 4000 6000]);

%%

clc;

compareTimes(tMin, tMax, ...
    xSolH, ...
    f, ...
    fg, ...
    fh, ...
    0, [], [], ...
    options, pSol,  xSigmaError, tSigmaError, [], [], [], ...
    {'ldl'} , ...
    [1000]);

%% 9
optionsCopy = struct(options);
optionsCopy.taskName = 'task_09';
if ( clearData )
    clc;
    close all;
end

r = 3.5;
m = 19;
pSol = [r; m];
optionsCopy.pSol = pSol;

pLb = [];
pUb = [];
%p = 1;
%pUb = Inf * ones(p, 1);
%pLb = - pUb;

x0 = 19.001;
p0 = [3.45; 19.1];

tau = 0.74;

tMin = 2 + tau;
tMax = 25;

delays = [ 0 tau ];

ddeOptions = ddeset('RelTol',1e-4,'AbsTol',1e-7,'InitialY', x0);
ddef = @(t, x, xLag, r, m) r * x * (1 - xLag / m);
ddeXHistory = 19;
ddeSol = dde23(ddef,tau,ddeXHistory,[tMin - tau, tMax],ddeOptions, pSol(1), pSol(2));

xSolH = @(t) interp1(ddeSol.x, ddeSol.y, t, 'spline', 'extrap');

xHistory = [];
[f, fg, fh] = ddeParamEst_checkUserInput('p_1 * x * ( 1 - x_1 / p_2 )', 2, 1, []);
% f = @(x, t, pSol) pSol(1) * x(1) * ( 1 - x(2) / pSol(2) );
% 
% fg = @(x, t, pSol) [ 
%         pSol(1) * ( 1 - x(2) / pSol(2) ), ...
%         - pSol(1) * x(1) / pSol(2), ...
%         x(1) * ( 1 + x(2) / pSol(2) ) , ...
%         pSol(1) * x(1) * x(2) / pSol(2) ^ 2 ];
% 
% fh = @(x, t, pSol) [
%     0,                                  1 - x(2) / pSol(2),            pSol(1) * x(2) / pSol(2) ^ 2; ...
%     1 - x(2) / pSol(2),                0,                              - x(1) * x(2) / pSol(2) ^ 2; ...
%     pSol(1) * x(2) / pSol(2) ^ 2,     - x(1) * x(2) / pSol(2) ^ 2    -2 * pSol(1) * x(1) * x(2)/ pSol(2) ^ 3];

ddeParamEst_TEST(...
    N,  tMin, tMax, ...
    xSolH, pSol, ...
    f, ...
    fg, ...
    fh, ...
    delays, xHistory, max(delays), ...
    optionsCopy, xSigmaError, tSigmaError, pLb, pUb, p0);

%%
compareTimes(tMin, tMax, ...
    xSolH, ...
    f, ...
    fg, ...
    fh, ...
    delays, xHistory, max(delays), ...
    options, pSol,  xSigmaError, tSigmaError, pLb, pUb, p0, ...
    {'ldl', 'bicg', 'fmincon'} , ...
    [100 500 1000]);

%% 10
options.taskName = 'task_10';

if ( clearData )
    clc;
    close all;
end

pSol = 3;

pUb = [];
pLb = [];
%p = 1;
%pUb = Inf * ones(p, 1);
%pLb = - pUb;

tau = 5;

tMin = 0;
tMax = 7;

delays = [ 0 tau ];

exampleoptOptions = ddeset('RelTol',1e-4,'AbsTol',1e-7);
exampleF = @(t,y,Z, pSol) pSol * Z;
exampleDelayF = @(t, pSol) 2;
exampleSolution = dde23(exampleF,tau,exampleDelayF,[tMin - tau, tMax],exampleoptOptions, pSol);

%xSolH = @(t) interp1q(exampleSolution.x', exampleSolution.y', t);
xSolH = @(t) interp1(exampleSolution.x, exampleSolution.y, t, 'spline', 'extrap');

% f = @(x, t, pSol) pSol * x(2);

xHistory = [];
    
% fg = @(x, t, pSol) [ 0, pSol, x(2) ];
% 
% fh = @(x, t, pSol) [
%     0,   1; ...
%     1,           0];

[f, fg, fh] = ddeParamEst_checkUserInput('p_1 * x_1', 1, 1, []);

ddeParamEst_TEST(...
    N, tMin, tMax, ...
    xSolH, pSol, ...
    f, ...
    fg, ...
    fh, ...
    delays, xHistory, max(delays), ...
    options, xSigmaError, tSigmaError, pLb, pUb);

%% PRE 11

r = 0.0257;
k = 0.566;
gamma = 1.623;
pSol = [r; k; gamma];

N0 = 1;
Nc = 5.2;

tau1 = 25;
tau2 = 30;
tau3 = 100;

tMin = 1970;
tMax = 2500;

ddeDelays = [ tau1, tau2, tau3 ];

ddeOptOptions = ddeset('RelTol',1e-7,'AbsTol',1e-10);

ddeExp = @(delays, p) exp( - p(2) * ( delays(3) - N0 ) );

K = @(delays, p) Nc + p(3) * (delays(2) - N0 ) * ddeExp(delays, p); 

ddef = @(t, x, delays, p) p(1) * (delays(1) ^ 2) * ( 1 - x / K(delays, p) );
ddeXHistory = @(t, pSol) (60100^2) * utils.arcctg((1995-t)/45) / (10^9);
ddeSol = dde23(ddef, ddeDelays, ddeXHistory,[tMin - max(ddeDelays), tMax], ddeOptOptions, pSol);

figure;
hold on;
grid on;
plot (ddeSol.x, ddeSol.y);

%% 11

%N = 50;

optionsCopy = struct(options);
optionsCopy.taskName = 'task_11';
optionsCopy.plotResult = true;
optionsCopy.plotExtResult = true;
%optionsCopy.checkHessian = true;
%optionsCopy.hessian_method = 'gauss-newton';
optionsCopy.hessian_method = 'newton';
optionsCopy.sqpOptions = struct(optionsCopy.sqpOptions);
optionsCopy.sqpOptions.stepMethodIterative = false;
optionsCopy.sqpOptions.stepMethod = 'ldl';
optionsCopy.sqpOptions.ldlsThreshold = 0.001;
%sqpOptions.algo_method        = 'quasi-Newton';
optionsCopy.sqpOptions.algo_method        = 'Newton';

optionsCopy.sqpOptions.algo_globalization = 'line-search';
%optionsCopy.sqpOptions.algo_globalization = 'unit step-size';

%optionsCopy.sqpOptions.stepMethod = 'bicg';
%optionsCopy.sqpOptions.iterativePrecondAlgorithmThresh = 0.0001;
%optionsCopy.sqpOptions.miter = 3;


optionsCopy.xTol = 10^-2;
optionsCopy.pTol = 10^-1;


if ( clearData )
    clc;
    close all;
end

r = 0.0257;
k = 0.566;
gamma = 1.623;
pSol = [r; k; gamma];
%p0 = 3 * [r; k; gamma];
p0 = 1.2 * [r; k; gamma];
optionsCopy.pSol = pSol;

N0 = 1;
Nc = 5.2;

pLb = [];
pUb = [];
%p = length(pSol);
%pLb = -2*p0;
%pUb = - pLb;

tau1 = 25;
tau2 = 30;
tau3 = 100;

tMin = 1970;
tMax = 2200;

optionsCopy.extTMax = tMax + 500;

ddeDelays = [ tau1, tau2, tau3 ];
delays = [ 0, tau1, tau2, tau3 ];

ddeOptOptions = ddeset('RelTol',1e-7,'AbsTol',1e-10);

ddeExp = @(delays, p) exp( - p(2) * ( delays(3) - N0 ) );

K = @(delays, p) Nc + p(3) * (delays(2) - N0 ) * ddeExp(delays, p); 

ddef = @(t, x, delays, p) p(1) * (delays(1) ^ 2) * ( 1 - x / K(delays, p) );

ddeXHistory = @(t, pSol) (62000^2) * utils.arcctg((2010-t)/45) / (10^9);

ddeSol = dde23(ddef, ddeDelays, ddeXHistory,[tMin - max(delays), tMax], ddeOptOptions, pSol);

%plot (ddeSol.x, ddeSol.y);

%xSolH = @(t) interp1q(ddeSolution.x', ddeSolution.y', t);
xSolH = @(t) interp1(ddeSol.x, ddeSol.y, t, 'spline', 'extrap');

pn = 3;
xHistory = xSolH;
userConstants = {{'Nc', Nc}, {'N0', N0}};
K = 'Nc + p_3 * (x_2 - N0) * exp( - p_2 * ( x_3 - N0 ) )';
[ f, fg, fh, fSym, fgSym fhSym] = ddeParamEst_checkUserInput(strcat('p_1 * (x_1 ^ 2) * (  1 - x / (', K, ')  )'), pn, 3, userConstants);

pretty(fSym);
pretty(fgSym);
pretty(fhSym);


xSigmaError = 0.0;
tSigmaError = 0.0;

ddeParamEst_TEST(...
    N,  tMin, tMax, ...
    xSolH, pSol,...
    f, ...
    fg, ...
    fh, ...
    delays, xHistory, max(delays), ...
    optionsCopy, xSigmaError, tSigmaError, pLb, pUb, p0);
%%
compareTimes(tMin, tMax, ...
    xSolH, ...
    f, ...
    fg, ...
    fh, ...
    delays, xHistory, max(delays), ...
    optionsCopy, pSol,  xSigmaError, tSigmaError, [], [], p0, ...
    {'ldl', 'fmincon', 'bicg'} , ...
    [100 150 200]);

%%

%{'fmincon', 'default', 'symrcm', 'amd', 'colamd', 'colperm', 'dmperm',
%'symamd'}

compareTimes(tMin, tMax, ...
    xSolH, ...
    f, ...
    fg, ...
    fh, ...
    delays, xHistory, max(delays), ...
    optionsCopy, pSol,  xSigmaError, tSigmaError, [], [], p0, ...
    {'fmincon', 'ldl', 'mldivide', 'bicg', 'cgs'} , ...
    [100 150 200 250 300 400]);

%% 12


% BROKEN DERIVATIVES

options.taskName = 'task_12';

if ( clearData )
    clc;
    close all;
end

r = 0.0257;
k = 0.566;
gamma = 1.623;
tau1 = 25;
tau2 = 30;
tau3 = 100;
ddeTheta = [r; k; gamma; tau1; tau2; tau3];

N0 = 1;
Nc = 5.2;

pLb = [];
pUb = [];%[0.2, 1, 3, 26, 31, 101]
%p = length(pSol);
%pLb = -Inf * ones(p, 1);
%pUb = - pLb;

pSol = [r; k; gamma; tau1; tau2; tau3];

tMin = 1500;
tMax = 2100;

ddeDelays = [ tau1, tau2, tau3 ];
delays = [ 0, -4, -5, -6 ];

ddeoptOptions = ddeset('RelTol',1e-7,'AbsTol',1e-10);

ddeExp = @(delays, pSol) exp( - pSol(2) * ( delays(3) - N0 ) );

K = @(delays, pSol) Nc + pSol(3) * (delays(2) - N0 ) * ddeExp(delays, pSol);

ddeFunction = @(t, x, delays, pSol) pSol(1) * (delays(1) ^ 2) * ( 1 - x / K(delays, pSol) );
ddeDelayFunction = @(t, pSol) 200 / ( 2025 - t) ;
ddeSolution = dde23(ddeFunction, ddeDelays, ddeDelayFunction,[tMin - max(getDelays(delays, pSol)), tMax], ddeoptOptions, ddeTheta);

%plot (ddeSolution.x, ddeSolution.y);

%xSolH = @(t) interp1q(ddeSolution.x', ddeSolution.y', t);
xSolH = @(t) interp1(ddeSolution.x, ddeSolution.y, t, 'spline', 'extrap');

f = @(x, t, pSol) ddeFunction(t, x(1), x(2:length(x)), pSol);

Kg_x = @(delays, pSol) ddeExp(delays, pSol) * pSol(3) * (1 - N0 + pSol(2) * delays(2));
Kg_theta1 = @(delays, pSol) 0;
Kg_theta2 = @(delays, pSol) - pSol(3) * ( delays(2) - N0 ) * ddeExp(delays, pSol) * (delays(3) - N0) ;
Kg_theta3 = @(delays, pSol) (delays(2) - N0) * ddeExp(delays, pSol);

fg_x = @(x, delays, pSol) delays(1) * ( 2 - ( ( ( 2 * x ) + delays(1) ) * K(delays, pSol) - delays(1) * x * Kg_x(delays, pSol) ) / K(delays, pSol) ^ 2 );
fg_theta1 = @(x, delays, pSol) delays(1) ^ 2 * ( 1 - ( x*K(delays, pSol) - pSol(1)*x*Kg_theta1(delays, pSol) ) / K(delays, pSol) ^ 2 );
fg_theta2 = @(x, delays, pSol) (pSol(1) * delays(1) ^ 2 * x * Kg_theta2(delays, pSol)) / K(delays, pSol) ^ 2;
fg_theta3 = @(x, delays, pSol) (pSol(1) * delays(1) ^ 2 * x * Kg_theta3(delays, pSol)) / K(delays, pSol) ^ 2;

fg0 = @(t, x, delays, pSol) [fg_x(x, delays, pSol), fg_theta1(x, delays, pSol), fg_theta2(x, delays, pSol), fg_theta3(x, delays, pSol)];
fg = @(x, t, pSol) fg0(t, x(1), x(2:length(x)), pSol);

xHistory = @(t) (62000^2) * utils.arcctg((2010-t)/45) / (10^9);
    
ddeParamEst_TEST(...
    N,  tMin, tMax, ...
    xSolH, ...
    f, ...
    fg, ...
    [], ...
    delays, xHistory, max(pUb(1, 4:length(pUb))), ...
    options, pSol, length(pSol), xSigmaError, tSigmaError, pLb, pUb);

%% 13


fileTaskName = 'population_china';
fileName = strcat('input/', fileTaskName, '.csv');

optionsCopy = struct(options);
optionsCopy.sqpOptions = struct(options.sqpOptions);
optionsCopy.taskName = strcat('task_13_', fileTaskName);
optionsCopy.xTol = 0.01;
optionsCopy.pTol = 0.1;
%optionsCopy.method = 'euler';
optionsCopy.method = 'backward-euler';
%optionsCopy.hessian_method = 'gauss-newton';
optionsCopy.hessian_method = 'newton';
%optionsCopy.maxNumberOfIterations = 7;
optionsCopy.sqp = true;
optionsCopy.plotExtResult = true;
optionsCopy.sqpOptions.stepMethodIterative = true;
%optionsCopy.sqpOptions.stepMethod = 'ldl';
optionsCopy.sqpOptions.stepMethod = 'bicg';

input = csvread(fileName, 1, 0);

[N, ~] = size(input);

optionsCopy.extTMax = input(end, 1) + 200;

from = find(input(1:end,1) == 1870);
tInput = input(from:N,1);
xInput = input(from:N,2) / 10 ^ 6;

if ( clearData )
    clc;
    close all;
end

r = 0.0257;
k = 0.566;
gamma = 1.623;
p0 = [r; k; gamma];

pSol = 10^2 * [0.001029007623221;
    0.346070601572724;
    1.525424216294689];
optionsCopy.pSol = pSol;


N0 = 0.23;
Nc = 1.2;

pLb = [];
pUb = [];
%p = length(pSol);
%pLb = -Inf * ones(p, 1);
%pUb = - pLb;

tau1 = 25;
tau2 = 30;
tau3 = 100;

delays = [ 0, tau1, tau2, tau3 ];

pn = length(p0);
userConstants = {{'Nc', Nc}, {'N0', N0}};
K = 'Nc + p_3 * (x_2 - N0) * exp( - p_2 * ( x_3 - N0 ) )';
[ f, fg, fh ] = ddeParamEst_checkUserInput(strcat('p_1 * (x_1 ^ 2) * (  1 - x / (', K, ')  )'), pn, 3, userConstants);

xHistory = [];%@(t) interp1(input(1:from,1), input(1:from,2) / 10^6, t, 'spline', 'extrap');%@(t) (62000^2) * utils.arcctg((2010-t)/45) / (10^9) / 6;
    
ddeParamEst(...
    tInput, xInput, ...
    f, ...
    fg, ...
    fh, ...
    delays, xHistory, max(delays), ...
    optionsCopy, pn, pLb, pUb, p0);

% saveas(h, strcat('output/', options.taskName, '_result'), 'png'); 

% save(strcat('output/', taskName, '.txt'), 'thetaResult', '-ascii');

%% 14

%fileTaskName = 'reduced2_population_japan';
%fileTaskName = 'reduced2_population_india';
%fileTaskName = 'reduced2_population_china_mod';
%fileTaskName = 'reduced2_population_china';
%fileTaskName = 'reduced_population_china';
fileTaskName = 'reduced_population_india';
fileName = strcat('input/', fileTaskName, '.csv');

optionsCopy = struct(options);
optionsCopy.sqpOptions = struct(options.sqpOptions);
optionsCopy.taskName = strcat('task_14_', fileTaskName);
optionsCopy.xTol = 0.001;
optionsCopy.pTol = 0.01;
optionsCopy.maxNumberOfIterations = 7;
optionsCopy.method = 'backward-euler';
optionsCopy.hessian_method = 'gauss-newton';
%optionsCopy.hessian_method = 'newton';
optionsCopy.sqp = true;
optionsCopy.plotExtResult = true;
optionsCopy.sqpOptions.stepMethodIterative = false;
optionsCopy.sqpOptions.stepMethod = 'ldls';
%optionsCopy.sqpOptions.stepMethod = 'bicgstab';
optionsCopy.sqpOptions.iterativePrecondAlgorithmThresh = 0.0001;


input = csvread(fileName, 1, 0);

[N, ~] = size(input);

m = 10^6;

optionsCopy.extTMax = input(end, 1) + 500;

from = find(input(1:end,1)==1870);
tInput = input(from:N,1);
xInput = input(from:N,2) / m;
% plot(tInput, xInput);

if ( clearData )
    clc;
    close all;
end

error = 0.5 * rand(3,1) - 0.5 / 2;

r = 0.243851567239265;
k = 1.096553138392189;
gamma = 0.765837675938692;
pSol = [r; k; gamma];
p0 = [r + r * error(1); k + r * error(2); gamma + r * error(3)];
p0 = [1; 1; 1];


N0 = 0.18 * (10 ^ 6) / m;
%N0 = 0.18 * 1.2 / 0.98 * (10 ^ 6) / m;
Nc = 0.98 * (10 ^ 6) / m;
%Nc = 1.2 * (10 ^ 6) / m;

pLb = [];
pUb = [];

tau1 = 25;
tau2 = 30;
tau3 = 100;

delays = [ 0, tau1, tau2, tau3 ];

pn = length(p0);
userConstants = {{'Nc', Nc}, {'N0', N0}};
K = 'Nc + p_3 * (x_2 - N0) * exp( - p_2 * ( x_3 - N0 ) )';
[ f, fg, fh ] = ddeParamEst_checkUserInput(strcat('p_1 * (x_1 ^ 2) * (  1 - x / (', K, ')  )'), pn, 3, userConstants);

[~, ~, ~] = ddeParamEst(...
    tInput, xInput, ...
    f, ...
    fg, ...
    fh, ...
    delays, [], max(delays), ...
    optionsCopy, pn, pLb, pUb, p0);

% saveas(h, strcat('output/', options.taskName, '_result'), 'png');
% 
% save(strcat('output/', options.taskName, '.txt'), 'thetaResult', '-ascii');

%% 14.1

%fileTaskName = 'reduced2_population_japan';
%fileTaskName = 'reduced2_population_india';
%fileTaskName = 'reduced2_population_china_mod';
fileTaskName = 'reduced_population_china';
fileName = strcat('input/', fileTaskName, '.csv');

optionsCopy = struct(options);
optionsCopy.sqpOptions = struct(options.sqpOptions);
optionsCopy.taskName = strcat('task-14.1-', fileTaskName);
optionsCopy.xTol = 0.0001;
optionsCopy.pTol = 0.01;
%optionsCopy.maxNumberOfIterations = 7;
optionsCopy.method = 'backward-euler';
optionsCopy.hessian_method = 'gauss-newton';
optionsCopy.sqp = true;
optionsCopy.plotExtResult = true;
optionsCopy.sqpOptions.stepMethodIterative = false;
optionsCopy.sqpOptions.stepMethod = 'ldls';
%optionsCopy.sqpOptions.stepMethod = 'bicg';

input = csvread(fileName, 1, 0);

[N, ~] = size(input);

m = 10^6;

optionsCopy.extTMax = input(end, 1) + 500;

from = find(input(1:end,1)==1871);
tInput = input(from:N,1);
xInput = input(from:N,2) / m;

if ( clearData )
    clc;
    close all;
end

r = 0.218687133409019;
k = 4.787953025884243;
gamma = 1.623;
pSol = [r; k];
p0 = 5 * [r; k];

N0 = 0.18 * (10 ^ 6) / m;
Nc = 0.98 * (10 ^ 6) / m;

pLb = [];
pUb = [];

tau1 = 25;
tau2 = 30;
tau3 = 100;

delays = [ 0, tau1, tau2, tau3 ];

pn = length(p0);
userConstants = {{'Nc', Nc}, {'N0', N0}, {'g', gamma}};
K = 'Nc + g * (x_2 - N0) * exp( - p_2 * ( x_3 - N0 ) )';
[ f, fg, fh ] = ddeParamEst_checkUserInput(strcat('p_1 * (x_1 ^ 2) * (  1 - x / (', K, ')  )'), pn, 3, userConstants);

[~, ~, ~] = ddeParamEst(...
    tInput, xInput, ...
    f, ...
    fg, ...
    fh, ...
    delays, [], max(delays), ...
    optionsCopy, pn, pLb, pUb, p0);

% saveas(h, strcat('output/', options.taskName, '_result'), 'png');
% 
% save(strcat('output/', options.taskName, '.txt'), 'thetaResult', '-ascii');

%% 14.2

fileTaskName = 'population_world';
fileName = strcat('input/', fileTaskName, '.csv');

optionsCopy = struct(options);
optionsCopy.sqpOptions = struct(options.sqpOptions);
optionsCopy.taskName = strcat('task-14.2-', fileTaskName);
optionsCopy.xTol = 0.001;
optionsCopy.pTol = 0.01;
%optionsCopy.maxNumberOfIterations = 7;
optionsCopy.method = 'backward-euler';
optionsCopy.hessian_method = 'gauss-newton';
%optionsCopy.hessian_method = 'newton';
optionsCopy.sqp = true;
optionsCopy.plotExtResult = true;
optionsCopy.sqpOptions.stepMethodIterative = false;
optionsCopy.sqpOptions.stepMethod = 'ldls';
%optionsCopy.sqpOptions.stepMethod = 'gmres';
%optionsCopy.sqpOptions.stepMethod = 'qmr';
%optionsCopy.sqpOptions.stepMethod = 'bicg';
%optionsCopy.sqpOptions.stepMethod = 'bicgstab';
optionsCopy.sqpOptions.iterativePrecondAlgorithmThresh = 0.01;

optionsCopy.sqpOptions.ldlsThreshold = 0.01;


input = csvread(fileName, 1, 0);

[N, ~] = size(input);

m = 10^6;

optionsCopy.extTMax = input(end, 1) + 200;

from = find(input(1:end,1)==1870);
tInput = input(from:N,1);
xInput = input(from:N,2) / m;
% plot(tInput, xInput);

if ( clearData )
    clc;
    close all;
end

error = zeros(3, 1);%0.5 * rand(3,1) - 0.5 / 2;

r = 0.018020578445818;
k = 2.256460994461019;
gamma = 4.574770144705490;
pSol = [r; k; gamma];
p0 = [r + r * error(1); k + k * error(2); gamma + gamma * error(3)];


N0 = 1 * (10 ^ 6) / m;
%N0 = 0.18 * 1.2 / 0.98 * (10 ^ 6) / m;
%{4; 4.5; 5; 6.5; 7; 10.5 }
Nc = 9 * (10 ^ 6) / m;
%Nc = 1.2 * (10 ^ 6) / m;

pLb = [];
pUb = [];

tau1 = 25;
tau2 = 30;
tau3 = 100;

delays = [ 0, tau1, tau2, tau3 ];

pn = length(p0);
userConstants = {{'Nc', Nc}, {'N0', N0}};
K = 'Nc + p_3 * (x_2 - N0) * exp( - p_2 * ( x_3 - N0 ) )';
[ f, fg, fh ] = ddeParamEst_checkUserInput(strcat('p_1 * (x_1 ^ 2) * (  1 - x / (', K, ')  )'), pn, 3, userConstants);

xHistory = [];%@(t, pSol) utils.populationKapitsa(t);

[~, ~, ~] = ddeParamEst(...
    tInput, xInput, ...
    f, ...
    fg, ...
    fh, ...
    delays, xHistory, max(delays), ...
    optionsCopy, pn, pLb, pUb, p0);

% saveas(h, strcat('output/', options.taskName, '_result'), 'png');
% 
% save(strcat('output/', options.taskName, '.txt'), 'thetaResult', '-ascii');
%%
optionsCopy.plotExtResult = false;
optionsCopy.plotResult = false;
for ldlThresh = [0.0001 0.001 0.01 0.02 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2]
    
    optionsCopy.sqpOptions.ldlsThreshold = ldlThresh;
    fprintf('### ldl thresh = %d', ldlThresh);
    
    [~, ~, ~] = ddeParamEst(...
        tInput, xInput, ...
        f, ...
        fg, ...
        fh, ...
        delays, xHistory, max(delays), ...
        optionsCopy, pn, pLb, pUb, p0);
    
end


%% 14.3

fileTaskName = 'population_world';
fileName = strcat('input/', fileTaskName, '.csv');

optionsCopy = struct(options);
optionsCopy.sqpOptions = struct(options.sqpOptions);
optionsCopy.taskName = strcat('task-14.2-', fileTaskName);
optionsCopy.xTol = 0.0001;
optionsCopy.pTol = 0.001;
%optionsCopy.maxNumberOfIterations = 7;
optionsCopy.method = 'backward-euler';
optionsCopy.hessian_method = 'gauss-newton';
%optionsCopy.hessian_method = 'newton';
optionsCopy.sqp = true;
optionsCopy.plotExtResult = true;
optionsCopy.sqpOptions.stepMethodIterative = false;
optionsCopy.sqpOptions.stepMethod = 'ldl';
%optionsCopy.sqpOptions.stepMethod = 'gmres';
%optionsCopy.sqpOptions.stepMethod = 'bicgstab';
%optionsCopy.sqpOptions.stepMethod = 'bicg';

input = csvread(fileName, 1, 0);

[N, ~] = size(input);

m = 10^6;

optionsCopy.extTMax = input(end, 1) + 500;

from = find(input(1:end,1)==1870);
tInput = input(from:N,1);
xInput = input(from:N,2) / m;
% plot(tInput, xInput);

if ( clearData )
    clc;
    close all;
end

error = zeros(3, 1);%0.5 * rand(3,1) - 0.5 / 2;

r = 0.018020578445818;
k = 2.256460994461019;
gamma = 1.623;
pSol = [r; k];
p0 = [r + r * error(1); k + k * error(2)];


N0 = 1 * (10 ^ 6) / m;
%N0 = 0.18 * 1.2 / 0.98 * (10 ^ 6) / m;
%{4; 4.5; 5; 6.5; 7; 10.5 }
Nc = 9 * (10 ^ 6) / m;
%Nc = 1.2 * (10 ^ 6) / m;

pLb = [];
pUb = [];

tau1 = 25;
tau2 = 30;
tau3 = 100;

delays = [ 0, tau1, tau2, tau3 ];

pn = length(p0);
userConstants = {{'Nc', Nc}, {'N0', N0}, {'g', gamma}};
K = 'Nc + g * (x_2 - N0) * exp( - p_2 * ( x_3 - N0 ) )';
[ f, fg, fh ] = ddeParamEst_checkUserInput(strcat('p_1 * (x_1 ^ 2) * (  1 - x / (', K, ')  )'), pn, 3, userConstants);

xHistory = [];%@(t, pSol) utils.populationKapitsa(t);

[~, ~, ~] = ddeParamEst(...
    tInput, xInput, ...
    f, ...
    fg, ...
    fh, ...
    delays, xHistory, max(delays), ...
    optionsCopy, pn, pLb, pUb, p0);

% saveas(h, strcat('output/', options.taskName, '_result'), 'png');
% 
% save(strcat('output/', options.taskName, '.txt'), 'thetaResult', '-ascii');
%%

if ( clearData )
    clc;
    close all;
end

r = 0.0257;
k = 0.566;
gamma = 1.623;
pSol = [r; k; gamma];

N0 = 0.18 * (10 ^ 6) / m;
Nc = 0.98 * (10 ^ 6) / m;

tau1 = 25;
tau2 = 30;
tau3 = 100;

delays = [ 0, tau1, tau2, tau3 ];
ndelays = length(delays);

pn = length(pSol);
userConstants = {{'Nc', Nc}, {'N0', N0}};
K = 'Nc + p_3 * (x_2 - N0) * exp( - p_2 * ( x_3 - N0 ) )';
f = ddeParamEst_checkUserInput(strcat('p_1 * (x_1 ^ 2) * (  1 - x / (', K, ')  )'), pn, ndelays, userConstants);
% p_1 * (x_1 ^ 2) * (  1 - x / (Nc + p_3 * (x_2 - 0.1) * exp( - p_2 * ( x_3 - 0.1 ) ))  )


tSpan = [1870, 2500];

% prepare DDE input data
ddeOptOptions = ddeset('AbsTol', min([options.xTol, options.pTol]));

% convert ddeParamEst input to dde23 input
ddeF = @(t, x, delays, p) f([x delays]', t, p);

ddeXHistory = @(t, pSol) utils.populationKapitsa(t);
% 3.844 * asin(1/sqrt(1+( (2001 - t) / 45)^2))

% do extrapolation
xdeResult = dde23(ddeF, delays(2:end), ddeXHistory, tSpan, ddeOptOptions, pSol);

t = linspace(tSpan(1), tSpan(2), 1000);
x = deval(xdeResult, t);

plot (t, x, '-g');
%% 15

fileTaskName = 'reduced2_population_china';
fileName = strcat('input/', fileTaskName, '.csv');

options.taskName = strcat('task_15_', fileTaskName);

input = csvread(fileName, 1, 0);

[N, ~] = size(input);

tInput = input(1:N,1);
xInput = input(1:N,2) / 10^6;

if ( clearData )
    clc;
    close all;
end


tau1 = 25;
tau2 = 30;
tau3 = 100;

N0 = 0.23;
Nc = 1.2;

delays = [ 0, -4, -5, -6 ];

pn = 6;
userConstants = {{'Nc', Nc}, {'N0', N0}};
K = 'Nc + p_3 * (x_2 - N0) * exp( - p_2 * ( x_3 - N0 ) )';
[ f, fg, ~ ] = ddeParamEst_checkUserInput(strcat('p_1 * (x_1 ^ 2) * (  1 - x / (', K, ')  )'), 3, 3, userConstants);

xHistory = [];
    
ddeParamEst(...
    tInput, xInput, ...
    f, ...
    fg, ...
    [], ...
    delays, xHistory, max([tau1, tau2, tau3]), ...
    options, 6);


end