function symBlockDecomposition_TEST( A, n )

tic;
[L, D] = symBlockDecomposition(A, n);
symBlockDecompositionTime = toc;
display(symBlockDecompositionTime);

tic;
[L_, D_, P_] = ldl(A);
ldlTime = toc;
display(ldlTime);

actual = L*D*L';
expected = A;

error = norm( expected - actual, inf );
display( error );
if ( error > 0.0000000001 ) 
    throw (MException ('AssertionError:ConditionFailed', 'Condition failed!'));
end

expected = P_'*A*P_;
actual = L_*D_*L_';
ldlError = norm( expected - actual, inf );
display( ldlError );   



end

