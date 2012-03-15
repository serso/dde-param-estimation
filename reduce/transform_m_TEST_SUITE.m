function transform_m_TEST_SUITE( )
%%TRANSFORM_M_TEST_SUITE Summary of this function goes here
%   Detailed explanation goes here

clc;

% Test 01
transform_m_TEST([1 2 3 4 5 6 7 8 9], [1 8 9], 1);
transform_m_TEST([1 2 3 4 5 6 7 8 9]', [1 8 9]', 1);
transform_m_TEST([1 2 3 4 5 6 7 8 9]', [1 7 8 9]', 2);
transform_m_TEST([1 2 3 4 5 6 7 8 9]', [1 6 7 8 9]', 3);
transform_m_TEST([1 2 3 4 5 6 7 8 9], [1 6 7 8 9], 3);

% Test 02
a = [   1 2 3 4; 
        5 6 7 8; 
        9 10 11 12;
        13 14 15 16];
expected = [1 3 4; 9 11 12; 13 15 16];

transform_m_TEST(a, expected, 1);

% Test 02
a = [   1 0 0 0 1 1 1;
        0 0 0 0 1 1 1;
        0 0 0 0 1 1 1;
        0 0 0 0 1 1 1;
        1 1 1 1 2 1 1;
        1 1 1 1 1 3 1;
        1 1 1 1 1 1 4;];
    
expected = [    1 1 1 1; 
                1 2 1 1; 
                1 1 3 1;
                1 1 1 4];

transform_m_TEST(a, expected, 2);

%%

end

