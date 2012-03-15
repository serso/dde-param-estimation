function getDelays_TEST( delays, theta, expected )
%%

actual = getDelays(delays, theta);

if ( actual ~= expected ) 
    throw (MException ('AssertionError:ConditionFailed', 'Expected result differes from actual!'));
end

%%
end

