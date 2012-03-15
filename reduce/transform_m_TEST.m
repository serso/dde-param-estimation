function transform_m_TEST( a, expected, p )
%TRANSFORM_M_TEST Summary of this function goes here
%   Detailed explanation goes here

actual = transform_m(a, p);

if ( norm(actual - expected, inf) > 0 ) 
    display(expected);
    display(actual);
    throw (MException ('AssertionError:ConditionFailed'));
end

end

