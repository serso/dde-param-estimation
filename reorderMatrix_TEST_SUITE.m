function reorderMatrix_TEST_SUITE()

%%
clc;

% M matrix with no delays (n ~ 10)

M = [2     0     0     0     0     0     0     0     0     0     0     0
    0     0     0     0     0     0     0     0     0     0     0     0
    0     0     2     0     0     0     0     0     0     0     0     0
    0     0     0     2     0     0     0     0     0     0     0     0
    0     0     0     0     2     0     0     0     0     0     0     0
    0     0     0     0     0     2     0     0     0     0     0     0
    0     0     0     0     0     0     2     0     0     0     0     0
    0     0     0     0     0     0     0     2     0     0     0     0
    0     0     0     0     0     0     0     0     2     0     0     0
    0     0     0     0     0     0     0     0     0     2     0     0
    0     0     0     0     0     0     0     0     0     0     2     0
    0     0     0     0     0     0     0     0     0     0     0     0];

ae = [   -1.0000    2.0000         0         0         0         0         0         0         0         0         0   -1.0000
    0   -1.0000    1.1111         0         0         0         0         0         0         0         0   -0.1111
    0         0   -1.0000    2.1111         0         0         0         0         0         0         0   -1.1111
    0         0         0   -1.0000    2.1111         0         0         0         0         0         0   -1.1111
    0         0         0         0   -1.0000    2.1111         0         0         0         0         0   -1.1111
    0         0         0         0         0   -1.0000    2.1111         0         0         0         0   -1.1111
    0         0         0         0         0         0   -1.0000    2.1111         0         0         0   -1.1111
    0         0         0         0         0         0         0   -1.0000    2.1111         0         0   -1.1111
    0         0         0         0         0         0         0         0   -1.0000    2.1111         0   -1.1111
    0         0         0         0         0         0         0         0         0   -1.0000    2.1111   -1.1111];

me = 10;

g = [       0
    0
    0.5516
    -0.0241
    0.0732
    0.0619
    -0.0811
    0.0377
    -0.0197
    -0.0187
    0
    0];

ce= [
    2.1967
    0.5412
    2.1700
    2.2096
    2.3359
    2.1565
    2.1610
    2.2792
    1.9879
    2.4147];

M0 = [M, ae'; ae, zeros(me)];
b0 = [-g; -ce];

M0 = sparse(M0);

reorderMatrix_TEST(M0, b0);

%% M matrix with no delays (n ~ 100)

load('reorderMatrix_TEST_01_M.mat', 'M0');
load('reorderMatrix_TEST_01_b.mat', 'b0');

reorderMatrix_TEST(M0, b0);

%% M matrix with delays (n ~ 10)

load('reorderMatrix_TEST_03_M.mat', 'M0');
load('reorderMatrix_TEST_03_b.mat', 'b0');

reorderMatrix_TEST(M0, b0);
%% M matrix with delays (n ~ 50)

load('reorderMatrix_TEST_02_M.mat', 'M0');
load('reorderMatrix_TEST_02_b.mat', 'b0');

reorderMatrix_TEST(M0, b0);


%%

end