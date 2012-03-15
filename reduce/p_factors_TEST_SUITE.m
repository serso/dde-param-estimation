function p_factors_TEST_SUITE()

%%  Test suite

clc;

%   Test 01

A = [   1 2 0 0 0 0 1;
        0 1 2 0 0 0 2;
        0 0 1 2 0 0 3;
        0 0 0 1 2 0 4;
        0 0 0 0 1 2 5];

b = [ 5; 4; 3; 2; 1 ];

p_factors_TEST(A, b, 1);

% Test 02

n = 10;
A = zeros(n - 1, n + 1);
b = zeros(n - 1, 1);
for i = 1 : n - 1
    A(i, i) = i;
    A(i, i + 1) = i + 1;
    b (i) = n - i;
end

p_factors_TEST(A, b, 1);

% Test 03

n = 10;
A = zeros(n - 1, n + 1);
b = zeros(n - 1, 1);
for i = 1 : n - 1
    A(i, i) = i;
    A(i, i + 1) = 1;
    
    A(i, n + 1) = 0.5;
     
    b (i) = 0.5;
end

p_factors_TEST(A, b, 1);

% Test 04

n = 10;

A = zeros(n - 1, n + 3);
b = zeros(n - 1, 1);
for i = 1 : n - 1
    A(i, i) = i;
    A(i, i + 1) = i + 1;
    
    A(i, n + 1) = 0.8;
    A(i, n + 2) = 0.7;
    A(i, n + 3) = 0.9;
    
    b (i) = 0.9;
end

p_factors_TEST(A, b, 3);


for n = [10 100 1000 2000 3000]
    display(n);
    
    A = zeros(n - 1, n + 3);
    b = zeros(n - 1, 1);
    for i = 1 : n - 1
        A(i, i) = i;
        A(i, i + 1) = i + 1;
        
        A(i, n + 1) = 0.8;
        A(i, n + 2) = 0.7;
        A(i, n + 3) = 0.9;
        
        b (i) = 0.9;
    end
    
    p_factors_TEST(A, b, 3);
end

% Test 05

for n = [10 100 1000 10000 20000]
    t = [0:1:n] / n;
    p = 1;

    addPath('../');

    A = getConstraintMatrix(t', n, p, [], [], 'euler');
    b = zeros(n - 1, 1);
    
    
    timerId = tic();
    
    p_factors_TEST(A, b, p);
    
    display(sprintf('Time: %0.5f s', toc(timerId)));
end

end

