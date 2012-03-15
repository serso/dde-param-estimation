function getDelays_TEST_SUITE()
%GETDELAYS_TEST_SUITE Test suite for getDelays method
%%

getDelays_TEST(0, [], 0);
getDelays_TEST(1, [], 1);
getDelays_TEST([23], [], [23]);
getDelays_TEST([0 2], [], [0 2]);
getDelays_TEST([0 2 25], [], [0 2 25]);

try
    getDelays_TEST([], [], []);
    throw (MException ('AssertionError:ConditionFailed', 'Expected result differes from actual!'));
catch ME
    % Get last segment of the error message identifier.
    idSegLast = regexp(ME.identifier, '(?<=:)\w+$', 'match');
    
    if strcmp(idSegLast, 'IllegalArgument')
        % ok, test succeded
    end
end

getDelays_TEST([0 2 25 -1], 3, [0 2 25 3]);
getDelays_TEST([0 2 25 -1 -2], [3 4], [0 2 25 3 4]);
getDelays_TEST([-1 2 25 -3 -2], [0 3 4], [0 2 25 4 3 ]);


%%
end

