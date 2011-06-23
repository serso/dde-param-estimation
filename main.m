function main()

%% intial clear data

clc;
close all;

% %% SQP
% 
% options = optimset('Algorithm', 'interior-point', 'MaxFunEvals', 3000);
% options = optimset(options, 'DerivativeCheck', 'off');
% options = optimset(options, 'FinDiffType', 'central');
% 
% %%
% grad = @(x) 2 * x;
% hess = @(x) 2;
% 
% thetaLb = 2;
% 
% [xResult, xStepsResult] = sqp( 1, grad, hess, [], [], 10, options, false, thetaLb, []);
% 
% display(xResult);
% display(xStepsResult);
% 
% %%
% 
% grad = @(x) [   2 * x(1) + 1 ; ...
%     2 * x(2) ];
% 
% hess = @(x) [   2,  0 ; ...
%     0,  2 ];
% 
% [xResult, xStepsResult] = sqp( 2, grad, hess, [], [], [10; 10], options, false, [], []);
% 
% display(xResult);
% display(xStepsResult);
% 
% %%
% 
% grad = @(x) [   2 * x(1) ; ...
%     2 * x(2) ];
% 
% hess = @(x) [   2,  0 ; ...
%     0,  2 ];
% 
% c = @(x)    x(1) + x(2);
% 
% jc = @(x)   [  1, 1 ];
% 
% [xResult, xStepsResult] = sqp( 2, grad, hess, c, jc, [10; 10], options, false, [], []);
% 
% display(xResult);
% display(xStepsResult);
% 
% %%
% 
% grad = @(x) [   2 * x(1) ; ...
%     2 * x(2) ];
% 
% hess = @(x) [   2,  0 ; ...
%     0,  2 ];
% 
% c = @(x)    1 / (x(1) + x(2)) - (x(1)^2 + x(2)^2);
% 
% jc = @(x)   [  - 1 / ( (x(1) + x(2)) * x(1) ) - 2 * x(1), - 1 / ( (x(1) + x(2)) * x(2) ) - 2 * x(2)];
% 
% [xResult, xStepsResult] = sqp( 2, grad, hess, c, jc, [10; 10], options, false, [], []);
% 
% display(xResult);
% display(xStepsResult);


%% Parameter estimation in ODE

debug = false;
showResult = true;
showIntermidiateResult = false;
clearData = true;

% methods = {'euler' 'backward_euler' 'box' 'rk4'};

% number of known points (i.e. values of function x(t))
N = 100;

xSigmaError = 0.05;
tSigmaError = 0.0;

options = optimset('Algorithm', 'interior-point', 'MaxFunEvals', 6000);
% options = optimset(options, 'SubproblemAlgorithm', 'cg');
options = optimset(options, 'DerivativeCheck', 'off');
options = optimset(options, 'FinDiffType', 'central');
% options = optimset(options, 'Algorithm', 'sqp');

method = 'backward_euler';

%% intial clear data

clc;
close all;
%% 1

taskName = 'task_01';    

if ( clearData )
    clc;
    close all;
end

tMin = 0;
tMax = 10;

theta = 2;

exactFun = @(t) theta - (theta - 1) * exp(-t);

fun = @(x, ~, theta) theta - x;

grad = @(~, ~, ~) [ -1, 1 ];
hess = @(x, t, theta) [ 0, 0; 0 0];

solveInManyIterations(...
    taskName, N,  tMin, tMax, ...
    exactFun, ...
    fun, ...
    grad, ...
    hess, ...
    0, [], [], ...
    method, options, theta, length(theta), xSigmaError, tSigmaError, debug, showResult, showIntermidiateResult);

%% 2

taskName = 'task_02';

if ( clearData )
    clc;
    close all;
end

tMin = 0;
tMax = 10;

theta = 2;

exactFun = @(t) theta - sin( t ) * exp(-t);

fun = @(x, t, theta) theta - x - exp( - t ) * cos (t);

grad = @(x, t, theta) [ -1, 1 ];
hess = @(x, t, theta) [ 0, 0; 0 0];

solveInManyIterations(...
    taskName, N,  tMin, tMax, ...
    exactFun, ...
    fun, ...
    grad, ...
    hess, ...
    0, [], [], ...
    method, options, theta, length(theta), xSigmaError, tSigmaError, debug, showResult, showIntermidiateResult);

%% 3

taskName = 'task_03';

if ( clearData )
    clc;
    close all;
end

tMin = 0.5;
tMax = 10;

theta = 2;

exactFun = @(t) theta * cos(t) / t;

fun = @(x, t, theta) -( theta * sin (t) / t + x / t  );

grad = @(x, t, theta) [ -1/t, -sin(t)/t ];
hess = @(x, t, theta) [ 0, 0; 0, 0];

solveInManyIterations(...
    taskName, N,  tMin, tMax, ...
    exactFun, ...
    fun, ...
    grad, ...
    hess, ...
    0, [], [], ...
    method, options, theta, length(theta), xSigmaError, tSigmaError, debug, showResult, showIntermidiateResult);

%% 4

taskName = 'task_04';

if ( clearData )
    clc;
    close all;
end

tMin = 0;
tMax = 10;

theta = [ 4; 2];

exactFun = @(t) theta(1) / theta(2) * ( 1 - exp( - theta(2) * t ));

fun = @(x, t, theta) theta(1) - theta(2) * x;

grad = @(x, t, theta) [ -theta(2), 1, -x ];
hess = @(x, t, theta) [
    0,  0, -1;
    0, 0, 0;
    -1, 0, 0];

solveInManyIterations(...
    taskName, N,  tMin, tMax, ...
    exactFun, ...
    fun, ...
    grad, ...
    hess, ...
    0, [], [], ...
    method, options, theta, length(theta), xSigmaError, tSigmaError, debug, showResult, showIntermidiateResult);

%% 5

taskName = 'task_05';

if ( clearData )
    clc;
    close all;
end

tMin = 0;
tMax = 10;

theta = [ 2; 1; 1; 1];

thetaLb = [ 0 0 0 0 ];
thetaUb = [ 8 5 5 5 ];
%p = 4;
%thetaUb = Inf * ones(p, 1);
%thetaLb = - thetaUb;

exactFun = @(t) theta(1) - theta(2) * sin( theta(3) * t ) * exp (- theta(4) * t);

fun = @(x, t, theta) theta(4) * (theta(1) - x) - theta(2) * theta(3) * cos(theta(3) * t ) * exp( - theta(4) * t );

grad = @(x, t, theta) [
    -theta(4)
    theta(4)
    -theta(3) * cos(theta(3) * t ) * exp( - theta(4) * t )
    -theta(2)  * exp( - theta(4) * t ) * (cos(theta(3) * t ) - theta(3) * sin(theta(3) * t) * t)
    (theta(1) - x) + t * theta(2) * theta(3) * cos(theta(3) * t ) * exp( - theta(4) * t )];

x_th4 = @(x, t, theta) -1;
th1_th4 = @(x, t, theta) 1;
th2_th3 = @(x, t, theta) exp( - theta(4) * t) * (theta(3) * t * sin(theta(3) * t) - cos(theta(3) * t) );
th2_th4 = @(x, t, theta) theta(3) * t * cos(theta(3) * t) * exp( - theta(4) * t);
th3_th3 = @(x, t, theta) theta(2) * exp( - theta(4) * t) * t * ( 2 * sin(theta(3) * t) + theta(3) * t * cos(theta(3) * t));
th3_th4 = @(x, t, theta) theta(2) * exp( - theta(4) * t) * t * ( cos(theta(3)) - theta(3) * t * sin(theta(3) * t));
th4_th4 = @(x, t, theta) theta(2) * theta(3) * t * t * cos (theta(3) * t) * exp( - theta(4) * t);

hess = @(x, t, theta) [
    0,                  0,                      0,                      0,                      x_th4(x, t, theta);
    0,                  0,                      0,                      0,                      th1_th4(x, t, theta);
    0,                  0,                      0,                      th2_th3(x, t, theta),   th2_th4(x, t, theta);
    0,                  0,                      th2_th3(x, t, theta),   th3_th3(x, t, theta),   th3_th4(x, t, theta);
    x_th4(x, t, theta), th1_th4(x, t, theta),   th2_th4(x, t, theta),   th3_th4(x, t, theta),   th4_th4(x, t, theta);];

solveInManyIterations(...
    taskName, N,  tMin, tMax, ...
    exactFun, ...
    fun, ...
    grad, ...
    hess, ...
    0, [], [], ...
    method, options, theta, length(theta), xSigmaError, tSigmaError, debug, showResult, showIntermidiateResult, thetaLb, thetaUb);

%% 6

taskName = 'task_06';

if ( clearData )
    clc;
    close all;
end

theta = 1;

tau = 2;

tMin = -tau;
tMax = 10;

%thetaLb = 0;
%thetaUb = 4;
p = 1;
thetaUb = Inf * ones(p, 1);
thetaLb = - thetaUb;

exactFun = @(t) cos( theta * pi * t / ( 2 * tau));

fun = @(x, t, theta) ( - theta * pi / (2 * tau)) * x(2);

delayF = []; %@(t) exactFun(theta, t);


grad = @(x, t, theta) [- theta * pi / (2 * tau), (- pi / (2 * tau)) * x(2)];

hess = @(x, t, theta) [
    0,                  - pi / (2 * tau); ...
    - pi / (2 * tau),     0];


delays = [ 0 tau ];

solveInManyIterations(...
    taskName, N,  tMin, tMax, ...
    exactFun, ...
    fun, ...
    grad, ...
    hess, ...
    delays, delayF, max(delays), ...
    method, options, theta, length(theta), xSigmaError, tSigmaError, debug, showResult, showIntermidiateResult, thetaLb, thetaUb);

%% 7
taskName = 'task_07';

if ( clearData )
    clc;
    close all;
end

theta = 1;
k = 1;
tau = pi / 2;
b = 1;

tMin = -tau;
tMax = 10;


%thetaLb = 0;
%thetaUb = 4;
p = 1;
thetaUb = Inf * ones(p, 1);
thetaLb = - thetaUb;

exactFun = @(t) b * sin ( k * theta * t ) - b * sin ( k * theta * tau ) / ( 1 - cos ( k * theta * tau )) * cos ( k * theta * t );

fun = @(x, t, theta) - theta * k * x(2);

delayF = []; %@(t) exactFun(theta, t);

grad = @(x, t, theta) [- theta * k, - k * x(2)];

hess = @(x, t, theta) [
    0,       - k; ...
    - k,     0];

delays = [ 0 tau ];

solveInManyIterations(...
    taskName, N,  tMin, tMax, ...
    exactFun, ...
    fun, ...
    grad, ...
    hess, ...
    delays, delayF, max(delays), ...
    method, options, theta, length(theta), xSigmaError, tSigmaError, debug, showResult, showIntermidiateResult, thetaLb, thetaUb);
%% 8
taskName = 'task_08';

if ( clearData )
    clc;
    close all;
end

theta = 1;

thetaUb = 2;
thetaLb = 0;
%p = 1;
%thetaUb = Inf * ones(p, 1);
%thetaLb = - thetaUb;

tau = 1;

tMin = tau;
tMax = 20;

delays = [ 0 tau ];

exampleF = @(t,y,Z,lambda) -lambda*Z*(1 + y);
exampleDelayF = @(t,lambda) t;
exampleSolution = dde23(exampleF,tau,exampleDelayF,[tMin - tau, tMax],[],theta);

exactFun = @(t) interp1(exampleSolution.x, exampleSolution.y, t, 'spline');
%exactFun = @(t) interp1q(exampleSolution.x', exampleSolution.y', t);

fun = @(x, t, theta) -theta * x(2) * (1 + x(1));

delayF = []; %@(t) t;

grad = @(x, t, theta) [-theta - theta * 2 * x(2), -x(2) * (1 + x(1))];

hess = @(x, t, theta) [
    - theta * 2,   -1 - 2 * x(2); ...
    -1 - 2 * x(2),           0];

solveInManyIterations(...
    taskName, N,  tMin, tMax, ...
    exactFun, ...
    fun, ...
    grad, ...
    hess, ...
    delays, delayF, max(delays), ...
    method, options, theta, length(theta), xSigmaError, tSigmaError, debug, showResult, showIntermidiateResult, thetaLb, thetaUb);

%% 9
taskName = 'task_09';

if ( clearData )
    clc;
    close all;
end

r = 3.5;
m = 19;
theta = [r; m];

thetaLb = [3 16];
thetaUb = [4 22];
%p = 1;
%thetaUb = Inf * ones(p, 1);
%thetaLb = - thetaUb;

x0 = 19.001;

tau = 0.74;

tMin = 2 + tau;
tMax = 25;

delays = [ 0 tau ];

exampleOptions = ddeset('RelTol',1e-4,'AbsTol',1e-7,'InitialY', x0);
exampleF = @(t,y,Z,r,m) r *y*(1 - Z/m);
exampleDelayF = 19;
exampleSolution = dde23(exampleF,tau,exampleDelayF,[tMin - tau, tMax],exampleOptions, theta(1), theta(2));

exactFun = @(t) interp1(exampleSolution.x, exampleSolution.y, t, 'spline', 'extrap');

fun = @(x, t, theta) theta(1) * x(1) * ( 1 - x(2) / theta(2) );

delayF = []; %@(t) 19;

grad = @(x, t, theta) [ 
        theta(1) * (1 - 2 * x(2)/ theta(2)), ...
        x(1) * ( 1 + x(2) / theta(2) ) , ...
        theta(1) * x(1) * x(2) / theta(2) ^ 2 ];

hess = @(x, t, theta) [
    0,                                  1 - x(2) / theta(2),            theta(1) * x(2) / theta(2) ^ 2; ...
    1 - x(2) / theta(2),                0,                              - x(1) * x(2) / theta(2) ^ 2; ...
    theta(1) * x(2) / theta(2) ^ 2,     - x(1) * x(2) / theta(2) ^ 2    -2 * theta(1) * x(1) * x(2)/ theta(2) ^ 3];

solveInManyIterations(...
    taskName, N,  tMin, tMax, ...
    exactFun, ...
    fun, ...
    [], ...
    [], ...
    delays, delayF, max(delays), ...
    method, options, theta, length(theta), xSigmaError, tSigmaError, debug, showResult, showIntermidiateResult, thetaLb, thetaUb, x0);

%% 10
taskName = 'task_10';

if ( clearData )
    clc;
    close all;
end

theta = 3;

thetaUb = 6;
thetaLb = 0;
%p = 1;
%thetaUb = Inf * ones(p, 1);
%thetaLb = - thetaUb;

tau = 5;

tMin = 0;
tMax = 7;

delays = [ 0 tau ];

exampleOptions = ddeset('RelTol',1e-4,'AbsTol',1e-7);
exampleF = @(t,y,Z, theta) theta * Z;
exampleDelayF = @(t, theta) 2;
exampleSolution = dde23(exampleF,tau,exampleDelayF,[tMin - tau, tMax],exampleOptions, theta);

%exactFun = @(t) interp1q(exampleSolution.x', exampleSolution.y', t);
exactFun = @(t) interp1(exampleSolution.x, exampleSolution.y, t, 'spline', 'extrap');

fun = @(x, t, theta) theta * x(2);

delayF = [];...@(t) t;
    
grad = @(x, t, theta) [ theta, x(2) ];

hess = @(x, t, theta) [
    0,   1; ...
    1,           0];

solveInManyIterations(...
    taskName, N, tMin, tMax, ...
    exactFun, ...
    fun, ...
    grad, ...
    hess, ...
    delays, delayF, max(delays), ...
    method, options, theta, length(theta), xSigmaError, tSigmaError, debug, showResult, showIntermidiateResult, thetaLb, thetaUb);

%% 11
taskName = 'task_11';

if ( clearData )
    clc;
    close all;
end

r = 0.0257;
k = 0.566;
gamma = 1.623;
theta = [r; k; gamma];

N0 = 1;
Nc = 5.2;

thetaLb = [0, 0, 0];
thetaUb = [0.2, 1, 3];
%p = length(theta);
%thetaLb = -Inf * ones(p, 1);
%thetaUb = - thetaLb;

tau1 = 25;
tau2 = 30;
tau3 = 100;

tMin = 0;
tMax = 500;

ddeDelays = [ tau1, tau2, tau3 ];
delays = [ 0, tau1, tau2, tau3 ];

ddeOptions = ddeset('RelTol',1e-7,'AbsTol',1e-10);

K = @(delays, theta) Nc + theta(3) * (delays(2) - N0 ) * exp( - theta(2) * ( delays(3) - N0 ) );

ddeFunction = @(t, x, delays, theta) theta(1) * (delays(1) ^ 2) * ( 1 - x / K(delays, theta) );
ddeDelayFunction = @(t, theta) 0.3;
ddeSolution = dde23(ddeFunction, ddeDelays, ddeDelayFunction,[tMin - max(delays), tMax], ddeOptions, theta);

% plot (ddeSolution.x, ddeSolution.y);

%exactFun = @(t) interp1q(ddeSolution.x', ddeSolution.y', t);
exactFun = @(t) interp1(ddeSolution.x, ddeSolution.y, t, 'spline', 'extrap');

fun = @(x, t, theta) ddeFunction(t, x(1), x(2:length(x)), theta);

delayF = [];...@(t) t;
    
solveInManyIterations(...
    taskName, N,  tMin, tMax, ...
    exactFun, ...
    fun, ...
    [], ...
    [], ...
    delays, delayF, max(delays), ...
    method, options, theta, length(theta), xSigmaError, tSigmaError, debug, showResult, showIntermidiateResult, thetaLb, thetaUb);

%% 12

taskName = 'task_12';

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

thetaLb = [0, 0, 0, 24, 29, 99];
thetaUb = [0.2, 1, 3, 26, 31, 101];
%p = length(theta);
%thetaLb = -Inf * ones(p, 1);
%thetaUb = - thetaLb;

theta = [r; k; gamma; tau1; tau2; tau3];

tMin = 1500;
tMax = 2100;

ddeDelays = [ tau1, tau2, tau3 ];
delays = [ 0, -4, -5, -6 ];

ddeOptions = ddeset('RelTol',1e-7,'AbsTol',1e-10);

K = @(delays, theta) Nc + theta(3) * (delays(2) - N0 ) * exp( - theta(2) * ( delays(3) - N0 ) );

ddeFunction = @(t, x, delays, theta) theta(1) * (delays(1) ^ 2) * ( 1 - x / K(delays, theta) );
ddeDelayFunction = @(t, theta) 200 / ( 2025 - t) ;
ddeSolution = dde23(ddeFunction, ddeDelays, ddeDelayFunction,[tMin - max(getDelays(delays, theta)), tMax], ddeOptions, ddeTheta);

%plot (ddeSolution.x, ddeSolution.y);

%exactFun = @(t) interp1q(ddeSolution.x', ddeSolution.y', t);
exactFun = @(t) interp1(ddeSolution.x, ddeSolution.y, t, 'spline', 'extrap');

fun = @(x, t, theta) ddeFunction(t, x(1), x(2:length(x)), theta);

delayF = [];...@(t) t;
    
solveInManyIterations(...
    taskName, N,  tMin, tMax, ...
    exactFun, ...
    fun, ...
    [], ...
    [], ...
    delays, delayF, max(thetaUb(1, 4:length(thetaUb))), ...
    method, options, theta, length(theta), xSigmaError, tSigmaError, debug, showResult, showIntermidiateResult, thetaLb, thetaUb);

%% 13


fileTaskName = 'population_china';
fileName = strcat('input/', fileTaskName, '.csv');

taskName = strcat('task_13_', fileTaskName);

input = csvread(fileName, 1, 0);

[N, ~] = size(input);

tInput = input(1:N,1);
xInput = input(1:N,2) / 10 ^ 6;

if ( clearData )
    clc;
    close all;
end

r = 0.0257;
k = 0.566;
gamma = 1.623;
theta0 = [r; k; gamma];

N0 = 0.23;
Nc = 1.2;

thetaLb = theta0 * 0;
thetaUb = theta0 * 1000;
%p = length(theta);
%thetaLb = -Inf * ones(p, 1);
%thetaUb = - thetaLb;

tau1 = 25;
tau2 = 30;
tau3 = 100;

delays = [ 0, tau1, tau2, tau3 ];

K = @(delays, theta) Nc + theta(3) * (delays(2) - N0 ) * exp( - theta(2) * ( delays(3) - N0 ) );

ddeFunction = @(t, x, delays, theta) theta(1) * (delays(1) ^ 2) * ( 1 - x / K(delays, theta) );

fun = @(x, t, theta) ddeFunction(t, x(1), x(2:length(x)), theta);

delayF = [];...@(t) t;
    
[~, ~, thetaResult]=solveInManyIterations1(...
    taskName, tInput, xInput, ...
    fun, ...
    [], ...
    [], ...
    delays, delayF, max(delays), ...
    method, options, length(theta0), debug, showResult, showIntermidiateResult, thetaLb, thetaUb, [], theta0);

h = plot(tInput, xInput, 'xr');
saveas(h, strcat('output/', taskName, '_result'), 'png'); 

save(strcat('output/', taskName, '.txt'), 'thetaResult', '-ascii');

%% 14

fileTaskName = 'reduced2_population_india';
fileName = strcat('input/', fileTaskName, '.csv');

taskName = strcat('task_14_', fileTaskName);

input = csvread(fileName, 1, 0);

[N, ~] = size(input);

tInput = input(1:N,1);
xInput = input(1:N,2) / 10^6;

if ( clearData )
    clc;
    close all;
end

r = 0.0257;
k = 0.566;
gamma = 1.623;
theta0 = [r; k; gamma];

N0 = 0.18;
Nc = 0.98;

thetaLb = theta0 * 0;
thetaUb = theta0 * 1000;
%p = length(theta);
%thetaLb = -Inf * ones(p, 1);
%thetaUb = - thetaLb;

tau1 = 25;
tau2 = 30;
tau3 = 100;

delays = [ 0, tau1, tau2, tau3 ];

K = @(delays, theta) Nc + theta(3) * (delays(2) - N0 ) * exp( - theta(2) * ( delays(3) - N0 ) );

ddeFunction = @(t, x, delays, theta) theta(1) * (delays(1) ^ 2) * ( 1 - x / K(delays, theta) );

fun = @(x, t, theta) ddeFunction(t, x(1), x(2:length(x)), theta);

delayF = [];...@(t) t;
    
[~, ~, thetaResult]=solveInManyIterations1(...
    taskName, tInput, xInput, ...
    fun, ...
    [], ...
    [], ...
    delays, delayF, max(delays), ...
    method, options, length(theta0), debug, showResult, showIntermidiateResult, thetaLb, thetaUb, [], theta0);

h = plot(tInput, xInput, 'xr');
saveas(h, strcat('output/', taskName, '_result'), 'png');

save(strcat('output/', taskName, '.txt'), 'thetaResult', '-ascii');

%% 15

fileTaskName = 'reduced2_population_china';
fileName = strcat('input/', fileTaskName, '.csv');

taskName = strcat('task_15_', fileTaskName);

input = csvread(fileName, 1, 0);

[N, ~] = size(input);

tInput = input(1:N,1);
xInput = input(1:N,2) / 10^6;

if ( clearData )
    clc;
    close all;
end

% theta = [r; k; gamma];


tau1 = 25;
tau2 = 30;
tau3 = 100;

N0 = 0.23;
Nc = 1.2;

thetaLb = [0, 0, 0, 24, 29, 99];
thetaUb = [0.2, 1, 3, 26, 31, 100];
%p = length(theta);
%thetaLb = -Inf * ones(p, 1);
%thetaUb = - thetaLb;

theta = [r; k; gamma; tau1; tau2; tau3];

delays = [ 0, -4, -5, -6 ];

K = @(delays, theta) Nc + theta(3) * (delays(2) - N0 ) * exp( - theta(2) * ( delays(3) - N0 ) );

ddeFunction = @(t, x, delays, theta) theta(1) * (delays(1) ^ 2) * ( 1 - x / K(delays, theta) );

fun = @(x, t, theta) ddeFunction(t, x(1), x(2:length(x)), theta);

delayF = [];...@(t) t;
    
solveInManyIterations1(...
    taskName, tInput, xInput, ...
    fun, ...
    [], ...
    [], ...
    delays, delayF, max(thetaUb(4:6 )), ...
    method, options, length(theta), debug, showResult, showIntermidiateResult, thetaLb, thetaUb);


end