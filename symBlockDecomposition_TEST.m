function symBlockDecomposition_TEST( A, n )

[L, D] = symBlockDecomposition(A, n);

actual = L*D*L';
expected = A;

error = norm( expected - actual, inf );
if ( error > 0.0000000001 ) 
    display( error );
    throw (MException ('AssertionError:ConditionFailed', 'Condition failed!'));
end


end

