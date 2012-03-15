function p_factors_TEST( A, b, p )
%TEST_P_FACTORS Function tests p_factors method

display('Testing...');

% display(A);
% display(b);

timeBefore = tic();
x = A \ b;
timeBefore = toc(timeBefore);
% display(x);
% display(A * x);

timeAfter = tic();
[P, P_add, new_A, new_b] = p_factors(A, b, p);

% display(P);
% display(P_add);

x_new = new_A \ new_b;
timeAfter = toc(timeAfter);

% display(x_new);
% display(P * x_new + P_add );

max_error = norm(P * x_new + P_add - x, inf);
display(max_error);

display(sprintf('Time before: %0.5f s', timeBefore));
display(sprintf('Time after: %0.5f s', timeAfter));

if max_error > 10^-10
     throw (MException ('AssertionError:ConditionFailed', 'Errors are too large - check the input data or algorithm!'));
end

end

