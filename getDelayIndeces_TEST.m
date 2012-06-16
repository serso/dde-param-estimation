function getDelayIndeces_TEST( t, delays, p, expected )
%%

actual = getDelayIndeces(t, delays, p);

if ( max(actual ~= expected) ) 
    display(actual);
    display(expected);
    throw (MException ('AssertionError:ConditionFailed', 'Expected result differes from actual!'));
end

%%
end

