function ddeParamEst_constraints_TEST_HESS()
%DDEPARAMEST_CONSTRAINTS_TEST_HESS Summary of this function goes here
%   Detailed explanation goes here

%%
addpath('derivset');

t = 1;

f = @(t) 2 - exp(-t);

display(hessian(f, t));

%%

end

