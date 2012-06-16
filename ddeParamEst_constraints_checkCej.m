function ddeParamEst_constraints_checkCej( cej, f, t, x, n, np, delays, xHistory, method, options )
%DDEPARAMEST_CONSTRAINTS_CHECKCEH Summary of this function goes here
%   Detailed explanation goes here
if ( isempty(f) )
    throw(MException('IllegalArgumentException:NotEmpty', 'f must be not empty!'));
end

    function result = fForFinDifCej (xParam)
        [~, result] = ddeParamEst_constraints(f, [], [], t, xParam, n, np, delays, xHistory, method, options, []);
    end

if ( n <= 110 )
    addpath('derivest');
    
    disp('Jacobian checking...');
    [finDifCej, ~] = jacobianest(@(xParam)fForFinDifCej(xParam), x);
    
    eps = 10E-10;
    finDifCej = utils.filterSmallElements(finDifCej, eps);
    
    if ( norm(finDifCej - cej, inf) > eps )
        %cej = finDifCej;
        figure;
        spy(cej);
        title('Analytically calculated jacobian');
        display(full(cej));
        
        figure;
        spy(finDifCej);
        title('Finite difference approximation of jacobian');
        display(finDifCej);
        
        jDiff = finDifCej - cej;
        jDiff = utils.filterSmallElements(jDiff, eps);
        figure;
        spy(jDiff);
        title('Difference between two jacobians');
        display(jDiff);
        
        throw(MException('AssertionError:ConditionFailed', 'Finite difference approximation of jacobian is different from analytically calculated!'));
    else
        disp('Jacobian check done!');
    end
else
    disp('Warning: Finite difference approximation could not be calculated due to big number of constraints: N must be less than 50');
end

end

