function getDelayIndeces_TEST_SUITE( )
%%

clc;

getDelayIndeces_TEST([1 2 3 4 5 6 7 8 9], [0 2 7], [], [1 3 8]);
getDelayIndeces_TEST([1 2 3 4 5 6 7 8 9], [0 2 7 -1], 4, [1 3 8 5]);
getDelayIndeces_TEST([1 2 3 4 5 6 7 8 9], 0, [], 1);


%%
end

