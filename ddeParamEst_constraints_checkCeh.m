function ddeParamEst_constraints_checkCeh( ceh, f, t, x, n, np, delays, xHistory, method, options, lm )
%DDEPARAMEST_CONSTRAINTS_CHECKCEH Summary of this function goes here
%   Detailed explanation goes here
if ( isempty(f) )
    throw(MException('IllegalArgumentException:NotEmpty', 'f must be not empty!'));
end

    function result = fForFinDifCeh (xParam)
        [~, ceRes] = ddeParamEst_constraints(f, [], [], t, xParam, n, np, delays, xHistory, method, options, lm);
        result = lm' * ceRes;
    end

if ( norm(lm, inf) > 0 )
    if ( n <= 60 )
        addpath('derivest');
        
        disp('Hessian checking...');
        [finDifCeh, ~] = hessian(@(xParam)fForFinDifCeh(xParam), x);
        
        eps = 10E-7;
        finDifCeh = utils.filterSmallElements(finDifCeh, eps);
        
        if ( norm(finDifCeh - ceh, inf) > eps )
            %ceh = finDifCeh;
            figure;
            spy(ceh);
            title('Analytically calculated hessian');
            display(full(ceh));
            
            figure;
            spy(finDifCeh);
            title('Finite difference approximation of hessian');
            display(finDifCeh);
            
            hDiff = finDifCeh - ceh;
            hDiff = utils.filterSmallElements(hDiff, eps);
            figure;
            spy(hDiff);
            title('Difference between two hessians');
            display(hDiff);
            
            throw(MException('AssertionError:ConditionFailed', 'Finite difference approximation of hessian is different from analytically calculated!'));
        else
            disp('Hessian check done!');
        end
    else
        disp('Warning: Finite difference approximation could not be calculated due to big number of constraints: N must be less than 50');
    end
else
    disp('lm == 0, hessian check skipped.');
end
    
end

