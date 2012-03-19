
function ddeParamEst_TEST_SUITE()

%% intial clear data

clc;
close all;

% %% SQP
% 
% optOptions = optimset('Algorithm', 'interior-point', 'MaxFunEvals', 3000);
% optOptions = optimset(optOptions, 'DerivativeCheck', 'off');
% optOptions = optimset(optOptions, 'FinDiffType', 'central');
% 
% %%
% fg = @(x) 2 * x;
% % fh = @(x) 2;
% 
% thetaLb = 2;
% 
% [xResult, xStepsResult] = sqp( 1, fg, fh, [], [], 10, optOptions, false, thetaLb, []);
% 
% display(xResult);
% display(xStepsResult);
% 
% %%
% 
% fg = @(x) [   2 * x(1) + 1 ; ...
%     2 * x(2) ];
% 
% fh = @(x) [   2,  0 ; ...
%     0,  2 ];
% 
% [xResult, xStepsResult] = sqp( 2, fg, fh, [], [], [10; 10], optOptions, false, [], []);
% 
% display(xResult);
% display(xStepsResult);
% 
% %%
% 
% fg = @(x) [   2 * x(1) ; ...
%     2 * x(2) ];
% 
% fh = @(x) [   2,  0 ; ...
%     0,  2 ];
% 
% c = @(x)    x(1) + x(2);
% 
% jc = @(x)   [  1, 1 ];
% 
% [xResult, xStepsResult] = sqp( 2, fg, fh, c, jc, [10; 10], optOptions, false, [], []);
% 
% display(xResult);
% display(xStepsResult);
% sqpoptOptions
% %%
% 
% fg = @(x) [   2 * x(1) ; ...
%     2 * x(2) ];
% 
% fh = @(x) [   2,  0 ; ...
%     0,  2 ];
% 
% c = @(x)    1 / (x(1) + x(2)) - (x(1)^2 + x(2)^2);
% 
% jc = @(x)   [  - 1 / ( (x(1) + x(2)) * x(1) ) - 2 * x(1), - 1 / ( (x(1) + x(2)) * x(2) ) - 2 * x(2)];
% 
% [xResult, xStepsResult] = sqp( 2, fg, fh, c, jc, [10; 10], optOptions, false, [], []);
% 
% display(xResult);
% display(xStepsResult);


%% Parameter estimation in ODE

showResult = true;
clearData = true;

% methods = {'euler' 'backward_euler' 'box' 'rk4'};

% number of known points (i.e. values of function x(t))
N = 100;

xSigmaError = 0.02;
tSigmaError = 0.0;

options.xTol = 10^-2;
options.thetaTol = 10^-2;

% optOptions = optimset('Algorithm', 'sqp');
optOptions = optimset('Algorithm', 'interior-point', 'MaxFunEvals', 6000);
% optOptions = optimset(optOptions, 'SubproblemAlgorithm', 'cg');
optOptions = optimset(optOptions, 'DerivativeCheck', 'off');
optOptions = optimset(optOptions, 'FinDiffType', 'central');
% optOptions = optimset(optOptions, 'Algorithm', 'sqp');

options.optOptions = optOptions;

options.sqp = true;
options.hessian_method = 'gauss-newton';
% options.hessian_method = 'newton';

% sqpOptions.algo_method        = 'quasi-Newton';
sqpOptions.algo_method        = 'Newton';
% sqpOptions.algo_globalization = 'line-search';
% sqpOptions.algo_globalization = 'unit step-size';
sqpOptions.stepMethod = 'symrcm';
sqpOptions.ldlsThreshold = 0.01;
sqpOptions.algo_globalization = 'line-search';
options.sqpOptions = sqpOptions;

%     sqpOptions.tol(1)  = tolopt(1);  % tolerance on the gradient of the Lagrangian
%     sqpOptions.tol(2)  = tolopt(2);  % tolerance on the feasibility
%     sqpOptions.tol(3)  = tolopt(3);  % tolerance on the complementarity

% options.sqpOptions = sqpOptions;

%method = 'backward_euler';
options.method = 'euler';

%options.method = method;
options.showResult = showResult; 

%% intial clear data

clc;
close all;
%% 1

options.taskName = 'task_01';    

if ( clearData )
    clc;    
    close all;
end

tMin = 0;
tMax = 10;

theta = 2;

x_sol = @(t) theta - (theta - 1) * exp(-t);

f = @(x, ~, theta) theta - x;

fg = @(~, ~, ~) [ -1, 1 ];
fh = @(x, t, theta) [ 0, 0; 0 0];

ddeParamEst_TEST(...
    N,  tMin, tMax, ...
    x_sol, ...
    f, ...
    fg, ...
    fh, ...
    0, [], [], ...
    options, theta, length(theta), xSigmaError, tSigmaError);

%%

compareTimes(tMin, tMax, ...
    x_sol, ...
    f, ...
    fg, ...
    fh, ...
    0, [], [], ...
    options, theta,  xSigmaError, tSigmaError, [], [], [], {'default', 'fmincon'});

%% 2

options.taskName = 'task_02';

if ( clearData )
    clc;
    close all;
end

tMin = 0;
tMax = 10;

theta = 2;      

x_sol = @(t) theta - sin( t ) * exp(-t);

f = @(x, t, theta) theta - x - exp( - t ) * cos (t);

fg = @(x, t, theta) [ -1, 1 ];
fh = @(x, t, theta) [ 0, 0; 0 0];

ddeParamEst_TEST(...
    N,  tMin, tMax, ...
    x_sol, ...
    f, ...
    fg, ...
    fh, ...
    0, [], [], ...
    options, theta, length(theta), xSigmaError, tSigmaError);

%% 3

options.taskName = 'task_03';

if ( clearData )
    clc;
    close all;
end

tMin = 0.5;
tMax = 10;

theta = 2;

x_sol = @(t) theta * cos(t) / t;

f = @(x, t, theta) -( theta * sin (t) / t + x / t  );

fg = @(x, t, theta) [ -1/t, -sin(t)/t ];
fh = @(x, t, theta) [ 0, 0; 0, 0];

ddeParamEst_TEST(...
    N,  tMin, tMax, ...
    x_sol, ...
    f, ...
    fg, ...
    fh, ...
    0, [], [], ...
    options, theta, length(theta), xSigmaError, tSigmaError);

%% 4

options.taskName = 'task_04';

if ( clearData )
    clc;
    close all;
end

tMin = 0;
tMax = 10;

theta = [ 4; 2];

x_sol = @(t) theta(1) / theta(2) * ( 1 - exp( - theta(2) * t ));

f = @(x, t, theta) theta(1) - theta(2) * x;

fg = @(x, t, theta) [ -theta(2), 1, -x ];
fh = @(x, t, theta) [
    0,  0, -1;
    0, 0, 0;
    -1, 0, 0];

ddeParamEst_TEST(...
    N,  tMin, tMax, ...
    x_sol, ...
    f, ...
    fg, ...
    fh, ...
    0, [], [], ...
    options, theta, length(theta), xSigmaError, tSigmaError);

%% 5

options.taskName = 'task_05';

if ( clearData )
    clc;
    close all;
end

tMin = 0;
tMax = 10;

theta = [ 2; 1; 1; 1];
theta0 = [ 2.1; 0.9; 0.9; 1.2];

thetaLb = [ 0 0 0 0 ];
thetaUb = [ 8 5 5 5 ];
%p = 4;
%thetaUb = Inf * ones(p, 1);
%thetaLb = - thetaUb;

x_sol = @(t) theta(1) - theta(2) * sin( theta(3) * t ) * exp (- theta(4) * t);

f = @(x, t, theta) theta(4) * (theta(1) - x) - theta(2) * theta(3) * cos(theta(3) * t ) * exp( - theta(4) * t );

fg = @(x, t, theta) [
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

fh = @(x, t, theta) [
    0,                  0,                      0,                      0,                      x_th4(x, t, theta);
    0,                  0,                      0,                      0,                      th1_th4(x, t, theta);
    0,                  0,                      0,                      th2_th3(x, t, theta),   th2_th4(x, t, theta);
    0,                  0,                      th2_th3(x, t, theta),   th3_th3(x, t, theta),   th3_th4(x, t, theta);
    x_th4(x, t, theta), th1_th4(x, t, theta),   th2_th4(x, t, theta),   th3_th4(x, t, theta),   th4_th4(x, t, theta);];

ddeParamEst_TEST(...
    N,  tMin, tMax, ...
    x_sol, ...
    f, ...
    fg, ...
    fh, ...
    0, [], [], ...
    options, theta, length(theta), xSigmaError, tSigmaError, thetaLb, thetaUb, theta0);

%% 6

N = 200;
options.showResult = true;
options.taskName = 'task_06';

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

x_sol = @(t) cos( theta * pi * t / ( 2 * tau));

f = @(x, t, theta) ( - theta * pi / (2 * tau)) * x(2);

delayF = []; %@(t) x_sol(theta, t);


fg = @(x, t, theta) [- theta * pi / (2 * tau), (- pi / (2 * tau)) * x(2)];

fh = @(x, t, theta) [
    0,                  - pi / (2 * tau); ...
    - pi / (2 * tau),     0];


delays = [ 0 tau ];

ddeParamEst_TEST(...
    N,  tMin, tMax, ...
    x_sol, ...
    f, ...
    fg, ...
    fh, ...
    delays, delayF, max(delays), ...
    options, theta, length(theta), xSigmaError, tSigmaError, thetaLb, thetaUb);
%%

compareTimes(tMin, tMax, ...
    x_sol, ...
    f, ...
    fg, ...
    fh, ...
    delays, delayF, max(delays), ...
    options, theta, xSigmaError, tSigmaError, thetaLb, thetaUb, 0.5, {'default', 'symrcm'}, [100 200 300 400 500]);

%% 7
options.taskName = 'task_07';

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

x_sol = @(t) b * sin ( k * theta * t ) - b * sin ( k * theta * tau ) / ( 1 - cos ( k * theta * tau )) * cos ( k * theta * t );

f = @(x, t, theta) - theta * k * x(2);

delayF = []; %@(t) x_sol(theta, t);

fg = @(x, t, theta) [- theta * k, - k * x(2)];

fh = @(x, t, theta) [
    0,       - k; ...
    - k,     0];

delays = [ 0 tau ];

ddeParamEst_TEST(...
    N,  tMin, tMax, ...
    x_sol, ...
    f, ...
    fg, ...
    fh, ...
    delays, delayF, max(delays), ...
    options, theta, length(theta), xSigmaError, tSigmaError, thetaLb, thetaUb);
%% 8
options.taskName = 'task_08';

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

x_sol = @(t) interp1(exampleSolution.x, exampleSolution.y, t, 'spline');
%x_sol = @(t) interp1q(exampleSolution.x', exampleSolution.y', t);

f = @(x, t, theta) -theta * x(2) * (1 + x(1));

delayF = []; %@(t) t;

fg = @(x, t, theta) [-theta - theta * 2 * x(2), -x(2) * (1 + x(1))];

fh = @(x, t, theta) [
    - theta * 2,   -1 - 2 * x(2); ...
    -1 - 2 * x(2),           0];

ddeParamEst_TEST(...
    N,  tMin, tMax, ...
    x_sol, ...
    f, ...
    fg, ...
    fh, ...
    delays, delayF, max(delays), ...
    options, theta, length(theta), xSigmaError, tSigmaError, thetaLb, thetaUb);
%%

compareTimes(tMin, tMax, ...
    x_sol, ...
    f, ...
    fg, ...
    fh, ...
    delays, delayF, max(delays), ...
    options, theta, xSigmaError, tSigmaError, thetaLb, thetaUb, 0.5, {'default', 'symrcm', 'fmincon'}, [100 200 300 400 500]);

%% 9
options.taskName = 'task_09';

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

exampleoptOptions = ddeset('RelTol',1e-4,'AbsTol',1e-7,'InitialY', x0);
exampleF = @(t,y,Z,r,m) r *y*(1 - Z/m);
exampleDelayF = 19;
exampleSolution = dde23(exampleF,tau,exampleDelayF,[tMin - tau, tMax],exampleoptOptions, theta(1), theta(2));

x_sol = @(t) interp1(exampleSolution.x, exampleSolution.y, t, 'spline', 'extrap');

f = @(x, t, theta) theta(1) * x(1) * ( 1 - x(2) / theta(2) );

delayF = []; %@(t) 19;

fg = @(x, t, theta) [ 
        theta(1) * (1 - 2 * x(2)/ theta(2)), ...
        x(1) * ( 1 + x(2) / theta(2) ) , ...
        theta(1) * x(1) * x(2) / theta(2) ^ 2 ];

fh = @(x, t, theta) [
    0,                                  1 - x(2) / theta(2),            theta(1) * x(2) / theta(2) ^ 2; ...
    1 - x(2) / theta(2),                0,                              - x(1) * x(2) / theta(2) ^ 2; ...
    theta(1) * x(2) / theta(2) ^ 2,     - x(1) * x(2) / theta(2) ^ 2    -2 * theta(1) * x(1) * x(2)/ theta(2) ^ 3];

ddeParamEst_TEST(...
    N,  tMin, tMax, ...
    x_sol, ...
    f, ...
    fg, ...
    fh, ...
    delays, delayF, max(delays), ...
    options, theta, length(theta), xSigmaError, tSigmaError, thetaLb, thetaUb, x0);

%% 10
options.taskName = 'task_10';

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

exampleoptOptions = ddeset('RelTol',1e-4,'AbsTol',1e-7);
exampleF = @(t,y,Z, theta) theta * Z;
exampleDelayF = @(t, theta) 2;
exampleSolution = dde23(exampleF,tau,exampleDelayF,[tMin - tau, tMax],exampleoptOptions, theta);

%x_sol = @(t) interp1q(exampleSolution.x', exampleSolution.y', t);
x_sol = @(t) interp1(exampleSolution.x, exampleSolution.y, t, 'spline', 'extrap');

f = @(x, t, theta) theta * x(2);

delayF = [];...@(t) t;
    
fg = @(x, t, theta) [ theta, x(2) ];

fh = @(x, t, theta) [
    0,   1; ...
    1,           0];

ddeParamEst_TEST(...
    N, tMin, tMax, ...
    x_sol, ...
    f, ...
    fg, ...
    fh, ...
    delays, delayF, max(delays), ...
    options, theta, length(theta), xSigmaError, tSigmaError, thetaLb, thetaUb);

%% 11

options.taskName = 'task_11';

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

ddeoptOptions = ddeset('RelTol',1e-7,'AbsTol',1e-10);

ddeExp = @(delays, theta) exp( - theta(2) * ( delays(3) - N0 ) );

K = @(delays, theta) Nc + theta(3) * (delays(2) - N0 ) * ddeExp(delays, theta); 

ddeFunction = @(t, x, delays, theta) theta(1) * (delays(1) ^ 2) * ( 1 - x / K(delays, theta) );
ddeDelayFunction = @(t, theta) 0.3;
ddeSolution = dde23(ddeFunction, ddeDelays, ddeDelayFunction,[tMin - max(delays), tMax], ddeoptOptions, theta);

% plot (ddeSolution.x, ddeSolution.y);

%x_sol = @(t) interp1q(ddeSolution.x', ddeSolution.y', t);
x_sol = @(t) interp1(ddeSolution.x, ddeSolution.y, t, 'spline', 'extrap');

f = @(x, t, theta) ddeFunction(t, x(1), x(2:length(x)), theta);

Kg_x = @(delays, theta) ddeExp(delays, theta) * theta(3) * (1 - N0 + theta(2) * delays(2));
Kg_theta1 = @(delays, theta) 0;
Kg_theta2 = @(delays, theta) - theta(3) * ( delays(2) - N0 ) * ddeExp(delays, theta) * (delays(3) - N0) ;
Kg_theta3 = @(delays, theta) (delays(2) - N0) * ddeExp(delays, theta);

fg_x = @(x, delays, theta) delays(1) * ( 2 - ( ( ( 2 * x ) + delays(1) ) * K(delays, theta) - delays(1) * x * Kg_x(delays, theta) ) / K(delays, theta) ^ 2 );
fg_theta1 = @(x, delays, theta) delays(1) ^ 2 * ( 1 - ( x*K(delays, theta) - theta(1)*x*Kg_theta1(delays, theta) ) / K(delays, theta) ^ 2 );
fg_theta2 = @(x, delays, theta) (theta(1) * delays(1) ^ 2 * x * Kg_theta2(delays, theta)) / K(delays, theta) ^ 2;
fg_theta3 = @(x, delays, theta) (theta(1) * delays(1) ^ 2 * x * Kg_theta3(delays, theta)) / K(delays, theta) ^ 2;

fg0 = @(t, x, delays, theta) [fg_x(x, delays, theta), fg_theta1(x, delays, theta), fg_theta2(x, delays, theta), fg_theta3(x, delays, theta)];
fg = @(x, t, theta) fg0(t, x(1), x(2:length(x)), theta);

delayF = [];...@(t) t;
    
ddeParamEst_TEST(...
    N,  tMin, tMax, ...
    x_sol, ...
    f, ...
    fg, ...
    [], ...
    delays, delayF, max(delays), ...
    options, theta, length(theta), xSigmaError, tSigmaError, thetaLb, thetaUb);

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

ddeoptOptions = ddeset('RelTol',1e-7,'AbsTol',1e-10);

ddeExp = @(delays, theta) exp( - theta(2) * ( delays(3) - N0 ) );

K = @(delays, theta) Nc + theta(3) * (delays(2) - N0 ) * ddeExp(delays, theta);

ddeFunction = @(t, x, delays, theta) theta(1) * (delays(1) ^ 2) * ( 1 - x / K(delays, theta) );
ddeDelayFunction = @(t, theta) 200 / ( 2025 - t) ;
ddeSolution = dde23(ddeFunction, ddeDelays, ddeDelayFunction,[tMin - max(getDelays(delays, theta)), tMax], ddeoptOptions, ddeTheta);

%plot (ddeSolution.x, ddeSolution.y);

%x_sol = @(t) interp1q(ddeSolution.x', ddeSolution.y', t);
x_sol = @(t) interp1(ddeSolution.x, ddeSolution.y, t, 'spline', 'extrap');

f = @(x, t, theta) ddeFunction(t, x(1), x(2:length(x)), theta);

Kg_x = @(delays, theta) ddeExp(delays, theta) * theta(3) * (1 - N0 + theta(2) * delays(2));
Kg_theta1 = @(delays, theta) 0;
Kg_theta2 = @(delays, theta) - theta(3) * ( delays(2) - N0 ) * ddeExp(delays, theta) * (delays(3) - N0) ;
Kg_theta3 = @(delays, theta) (delays(2) - N0) * ddeExp(delays, theta);

fg_x = @(x, delays, theta) delays(1) * ( 2 - ( ( ( 2 * x ) + delays(1) ) * K(delays, theta) - delays(1) * x * Kg_x(delays, theta) ) / K(delays, theta) ^ 2 );
fg_theta1 = @(x, delays, theta) delays(1) ^ 2 * ( 1 - ( x*K(delays, theta) - theta(1)*x*Kg_theta1(delays, theta) ) / K(delays, theta) ^ 2 );
fg_theta2 = @(x, delays, theta) (theta(1) * delays(1) ^ 2 * x * Kg_theta2(delays, theta)) / K(delays, theta) ^ 2;
fg_theta3 = @(x, delays, theta) (theta(1) * delays(1) ^ 2 * x * Kg_theta3(delays, theta)) / K(delays, theta) ^ 2;

fg0 = @(t, x, delays, theta) [fg_x(x, delays, theta), fg_theta1(x, delays, theta), fg_theta2(x, delays, theta), fg_theta3(x, delays, theta)];
fg = @(x, t, theta) fg0(t, x(1), x(2:length(x)), theta);

delayF = [];...@(t) t;
    
ddeParamEst_TEST(...
    N,  tMin, tMax, ...
    x_sol, ...
    f, ...
    fg, ...
    [], ...
    delays, delayF, max(thetaUb(1, 4:length(thetaUb))), ...
    options, theta, length(theta), xSigmaError, tSigmaError, thetaLb, thetaUb);

%% 13


fileTaskName = 'population_china';
fileName = strcat('input/', fileTaskName, '.csv');

options.taskName = strcat('task_13_', fileTaskName);

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


ddeExp = @(delays, theta) exp( - theta(2) * ( delays(3) - N0 ) );

delays = [ 0, tau1, tau2, tau3 ];

N1 = @(delays, theta) theta(3) * (delays(2) - N0 );
N2 = @(delays, theta) exp( - theta(2) * ( delays(3) - N0 ) );
K = @(delays, theta) Nc + N1(delays, theta) * N2(delays, theta);

ddeFunction = @(t, x, delays, theta) theta(1) * (delays(1) ^ 2) * ( 1 - x / K(delays, theta) );

f = @(x, t, theta) ddeFunction(t, x(1), x(2:length(x)), theta);

Kg_x = @(delays, theta) N2(delays, theta) * theta(3) * (1 - theta(2) * (delays(2) - N0) );
% Kg_theta1 = @(delays, theta) 0;
Kg_theta2 = @(delays, theta) - N1(delays, theta) * N2(delays, theta) * (delays(2) - N0) ;
Kg_theta3 = @(delays, theta) (delays(2) - N0) * N2(delays, theta);

fg_x = @(x, delays, theta) delays(1) * ( 2 - ( ( ( 2 * x ) + delays(1) ) * K(delays, theta) - delays(1) * x * Kg_x(delays, theta) ) / K(delays, theta) ^ 2 );
fg_theta1 = @(x, delays, theta) delays(1) ^ 2 * ( 1 - x/K(delays, theta));
fg_theta2 = @(x, delays, theta) (theta(1) * delays(1) ^ 2 * x * Kg_theta2(delays, theta)) / K(delays, theta) ^ 2;
fg_theta3 = @(x, delays, theta) (theta(1) * delays(1) ^ 2 * x * Kg_theta3(delays, theta)) / K(delays, theta) ^ 2;

fg0 = @(t, x, delays, theta) [fg_x(x, delays, theta), fg_theta1(x, delays, theta), fg_theta2(x, delays, theta), fg_theta3(x, delays, theta)];
fg = @(x, t, theta) fg0(t, x(1), x(2:length(x)), theta);

delayF = [];...@(t) t;
    
[~, ~, thetaResult]=ddeParamEst(...
    tInput, xInput, ...
    f, ...
    fg, ...
    [], ...
    delays, delayF, max(delays), ...
    options, length(theta0), thetaLb, thetaUb, theta0, []);

h = plot(tInput, xInput, 'xr');
saveas(h, strcat('output/', options.taskName, '_result'), 'png'); 

save(strcat('output/', taskName, '.txt'), 'thetaResult', '-ascii');

%% 14

fileTaskName = 'reduced2_population_china';
% fileTaskName = 'reduced2_population_india';
fileName = strcat('input/', fileTaskName, '.csv');

options.taskName = strcat('task_14_', fileTaskName);

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

thetaLb = theta0 - theta0 * 0.5;
thetaUb = theta0 + theta0 * 10;
%p = length(theta);
%thetaLb = -Inf * ones(p, 1);
%thetaUb = - thetaLb;

tau1 = 25;
tau2 = 30;
tau3 = 100;

delays = [ 0, tau1, tau2, tau3 ];


N1 = @(delays, theta) theta(3) * (delays(2) - N0 );
N2 = @(delays, theta) exp( - theta(2) * ( delays(3) - N0 ) );
K = @(delays, theta) Nc + N1(delays, theta) * N2(delays, theta);

ddeFunction = @(t, x, delays, theta) theta(1) * (delays(1) ^ 2) * ( 1 - x / K(delays, theta) );

f = @(x, t, theta) ddeFunction(t, x(1), x(2:length(x)), theta);

Kg_x = @(delays, theta) N2(delays, theta) * theta(3) * (1 - theta(2) * (delays(2) - N0) );
% Kg_theta1 = @(delays, theta) 0;
Kg_theta2 = @(delays, theta) - N1(delays, theta) * N2(delays, theta) * (delays(2) - N0) ;
Kg_theta3 = @(delays, theta) (delays(2) - N0) * N2(delays, theta);

fg_x = @(x, delays, theta) delays(1) * ( 2 - ( ( ( 2 * x ) + delays(1) ) * K(delays, theta) - delays(1) * x * Kg_x(delays, theta) ) / K(delays, theta) ^ 2 );
fg_theta1 = @(x, delays, theta) delays(1) ^ 2 * ( 1 - x/K(delays, theta));
fg_theta2 = @(x, delays, theta) (theta(1) * delays(1) ^ 2 * x * Kg_theta2(delays, theta)) / K(delays, theta) ^ 2;
fg_theta3 = @(x, delays, theta) (theta(1) * delays(1) ^ 2 * x * Kg_theta3(delays, theta)) / K(delays, theta) ^ 2;

fg0 = @(t, x, delays, theta) [fg_x(x, delays, theta), fg_theta1(x, delays, theta), fg_theta2(x, delays, theta), fg_theta3(x, delays, theta)];
fg = @(x, t, theta) fg0(t, x(1), x(2:length(x)), theta);

delayF = [];...@(t) t;
    
[~, ~, thetaResult]=ddeParamEst(...
    tInput, xInput, ...
    f, ...
    fg, ...
    [], ...
    delays, delayF, max(delays), ...
    options, length(theta0), thetaLb, thetaUb, theta0, []);

h = plot(tInput, xInput, 'xr');
saveas(h, strcat('output/', options.taskName, '_result'), 'png');

save(strcat('output/', options.taskName, '.txt'), 'thetaResult', '-ascii');

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

theta0 = [r; k; gamma];


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
N1 = @(delays, theta) theta(3) * (delays(2) - N0 );
N2 = @(delays, theta) exp( - theta(2) * ( delays(3) - N0 ) );
K = @(delays, theta) Nc + N1(delays, theta) * N2(delays, theta);

ddeFunction = @(t, x, delays, theta) theta(1) * (delays(1) ^ 2) * ( 1 - x / K(delays, theta) );

f = @(x, t, theta) ddeFunction(t, x(1), x(2:length(x)), theta);

Kg_x = @(delays, theta) N2(delays, theta) * theta(3) * (1 - theta(2) * (delays(2) - N0) );
% Kg_theta1 = @(delays, theta) 0;
Kg_theta2 = @(delays, theta) - N1(delays, theta) * N2(delays, theta) * (delays(2) - N0) ;
Kg_theta3 = @(delays, theta) (delays(2) - N0) * N2(delays, theta);

fg_x = @(x, delays, theta) delays(1) * ( 2 - ( ( ( 2 * x ) + delays(1) ) * K(delays, theta) - delays(1) * x * Kg_x(delays, theta) ) / K(delays, theta) ^ 2 );
fg_theta1 = @(x, delays, theta) delays(1) ^ 2 * ( 1 - x/K(delays, theta));
fg_theta2 = @(x, delays, theta) (theta(1) * delays(1) ^ 2 * x * Kg_theta2(delays, theta)) / K(delays, theta) ^ 2;
fg_theta3 = @(x, delays, theta) (theta(1) * delays(1) ^ 2 * x * Kg_theta3(delays, theta)) / K(delays, theta) ^ 2;

fg0 = @(t, x, delays, theta) [fg_x(x, delays, theta), fg_theta1(x, delays, theta), fg_theta2(x, delays, theta), fg_theta3(x, delays, theta)];
fg = @(x, t, theta) fg0(t, x(1), x(2:length(x)), theta);

delayF = [];...@(t) t;
    
ddeParamEst(...
    tInput, xInput, ...
    f, ...
    fg, ...
    [], ...
    delays, delayF, max(thetaUb(4:6 )), ...
    options, length(theta), thetaLb, thetaUb, theta0, []);


end