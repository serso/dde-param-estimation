function getDelayIndeces_TEST_SUITE( )
%%

clc;

getDelayIndeces_TEST([1 2 3 4 5 6 7 8 9], [0 2 7], [], [1 3 8]);
getDelayIndeces_TEST([1 2 3 4 5 6 7 8 9], [0 3.2 7], [], [1 4 8]);
getDelayIndeces_TEST([1 2 3 4 5 6 7 8 9], [0 3.51 7], [], [1 5 8]);
getDelayIndeces_TEST([1 2 3 4 5 6 7 8 9], [0 3 7], [], [1 4 8]);
getDelayIndeces_TEST([1 2 3 4 5 6 7 8 9], [0 2 7 -1], 4, [1 3 8 5]);
getDelayIndeces_TEST([1 2 3 4 5 6 7 8 9], 0, [], 1);
getDelayIndeces_TEST(linspace(0, 4.5, 10), [0 2], [], [1 5]);
getDelayIndeces_TEST(linspace(0, 10, 26), [0 6], [], [1 16]);


%%
end

