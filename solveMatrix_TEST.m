function solveMatrix_TEST( )
%%SOLVEMATRIX_TEST Summary of this function goes here
%   Detailed explanation goes here
clc;

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
%%
[L,D,P] = ldl(M0);
% x=(P*S')*(L'\(D\(L\((S*P)*b0))));
x= P * (L'\(D\(L\(P' * b0))));

solError = norm(M0*x-b0, inf);
display(solError);
if  solError > 10^-10
    throw (MException ('AssertionError:ConditionFailed', 'Not a solution!'));
end

%%
[L,D,P,S] = ldl(M0);
% x=(P*S')*(L'\(D\(L\((S*P)*b0))));
spy(S);
x= (S*P) * (L'\(D\(L\( (P'*S) * b0))));

solError = norm(M0*x-b0, inf);
display(solError);
if  solError > 10^-10
    throw (MException ('AssertionError:ConditionFailed', 'Not a solution!'));
end

%%
[L, U, P, Q] = lu(M0);

y = L \ (P * b0);
x = Q * (U \ y);

solError = norm(M0*x-b0, inf);
if  solError > 10^-10
    display(solError);
    throw (MException ('AssertionError:ConditionFailed', 'Not a solution!'));
end

%%
p = symrcm(M0);
M0p = M0(p,p);
spy(M0, 'b');
hold on;
spy(M0p, 'r');
b0p = b0(p);

xp = M0p \ b0p;
[~, rp] = sort(p);
x = xp(rp);

solError = norm(M0*x-b0, inf);
if  solError > 10^-10
    display(solError);
    throw (MException ('AssertionError:ConditionFailed', 'Not a solution!'));
end

%%
p = amd(M0);
M0p = M0(p,p);
spy(M0, 'b');
hold on;
spy(M0p, 'r');
b0p = b0(p);

xp = M0p \ b0p;
[~, rp] = sort(p);
x = xp(rp);

solError = norm(M0*x-b0, inf);
if  solError > 10^-10
    display(solError);
    throw (MException ('AssertionError:ConditionFailed', 'Not a solution!'));
end

%%
p = colamd(M0);
M0p = M0(p,p);
spy(M0, 'b');
hold on;
spy(M0p, 'r');
b0p = b0(p);

xp = M0p \ b0p;
[~, rp] = sort(p);
x = xp(rp);

solError = norm(M0*x-b0, inf);
if  solError > 10^-10
    display(solError);
    throw (MException ('AssertionError:ConditionFailed', 'Not a solution!'));
end

%%
p = colperm(M0);
M0p = M0(p,p);
spy(M0, 'b');
hold on;
spy(M0p, 'r');
b0p = b0(p);

xp = M0p \ b0p;
[~, rp] = sort(p);
x = xp(rp);

solError = norm(M0*x-b0, inf);
if  solError > 10^-10
    display(solError);
    throw (MException ('AssertionError:ConditionFailed', 'Not a solution!'));
end

%%
p = dmperm(M0);
M0p = M0(p,p);
spy(M0, 'b');
hold on;
spy(M0p, 'r');
b0p = b0(p);

xp = M0p \ b0p;
[~, rp] = sort(p);
x = xp(rp);

solError = norm(M0*x-b0, inf);
if  solError > 10^-10
    display(solError);
    throw (MException ('AssertionError:ConditionFailed', 'Not a solution!'));
end

%%

p = symamd(M0);
M0p = M0(p,p);
spy(M0, 'b');
hold on;
spy(M0p, 'r');
b0p = b0(p);

xp = M0p \ b0p;
[~, rp] = sort(p);
x = xp(rp);

solError = norm(M0*x-b0, inf);
if  solError > 10^-10
    display(solError);
    throw (MException ('AssertionError:ConditionFailed', 'Not a solution!'));
end

%%
end

