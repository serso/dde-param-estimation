function getDelayIndeces_TEST( t, delays, theta, expected )
%%

actual = getDelayIndeces(t, delays, theta);

if ( actual ~= expected ) 
    display(actual);
    display(expected);
    throw (MException ('AssertionError:ConditionFailed', 'Expected result differes from actual!'));
end

%%
end

