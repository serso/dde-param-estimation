function  [x,lm,info] = sqp_reduce( f, ce, cs_jacobian, x0 )
%%SQP_REDUCE Summary of this function goes here
%   Detailed explanation goes here


addpath(fullfile(pwd, '../sqp/src'));

% default options

% options.algo_method        = 'quasi-Newton';
options.algo_method        = 'Newton';

options.algo_globalization = 'line-search';
% options.algo_globalization = 'unit step-size';

    function [outdic, f_res, ci_res, ce_res, cs_res, g_res, ai_res, ae_res ] = mysimul ( indic, x, lm )
       
        f_res = [];
        ci_res = [];
        ce_res = [];
        cs_res = [];
        g_res = [];
        ai_res = [];
        ae_res = [];
        
        display(indic);
        
        if indic <= 1
            outdic = 0;
        elseif indic <= 2
            outdic = 0;
            f_res = f(x);
            ce_res = ce(x);
        elseif indic == 5
            outdic = 0;
            f_res = cs_jacobian(x);
        else 
            outdic = -2;
        end
    end

[n, ~] = size(x0);
lb = -inf * ones(n, 1);
ub = inf * ones(n, 1);

[x,lm,info] = sqplab(@mysimul, x0, [], lb, ub, options);

%%
end

