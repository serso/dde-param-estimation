function [ P, P_add, A_, b_] = p_factors( A, b, p, method )

%% REDUCE
% Function returns factors P, P_add which reduces the number of variables
% of equations A * x = b
% where A is a matrix of special form:
%
% A =
%     a1     b1     0      0      ...     0      g1_1       g2_1        gp_1
%     0      a2     b2     0      ...     0      g1_2       g2_2        gp_2
%     0      0      a3     b3     ...     0      g1_3       g2_3        gp_3
%     0      0      0      a4     ...     0      g1_4       g2_4        gp_4
%     0      0      0      0      ...     0      g1_5       g2_5        gp_5
%     ...    ...    ...    ...    ...     ...    ...        ...         ...
%     0      0      0      0      0       b(n-1) g1_(n-1)   g2_(n_1)    gp_(n-1)
%
%  All b(i) must be not zeros.
%
% ] d(i) = - a(i)/b(i)
%   e(i, j) = - g(i, j) / b(i)
%   f(i) = b_(i) / b(i) NOTE: b_(i) is an i-th element of input vector b 
%
% Next recurrent formalae were retreived:
%
% x(i+1) = x(1) * D(i+1) + x(n+1) * E(i+1) + F(i+1)
%
% where
%
% D(i+1) = ∏( d(j), j, 1, i ) = d(i) * D(i) * P
% E(i+1) = ∑( e(j) * ∏(d(k), k, j+1, i), j, 1, i) = e(i) + d(i) * E(i)
% D(i+1) = ∑( f(j) * ∏(d(k), k, j+1, i), j, 1, i) = f(i) + d(i) * F(i)

if ( nargin < 4 )
    method = 'fast';
end

[rows, ~] = size(A);
n = rows + 1;

d = zeros(n, 1);
e = zeros(n, p);
f = zeros(n, 1);

% x = P * x_ + P_add
P = zeros (n + p, 3);

P(1, 1:3) = [1 0 0];
P(n, 1:3) = [0 1 0];

P_add = zeros (n + p, 1);

P_add(1) = 0;
P_add(n) = 0;

for i = 1:p
    P(n + i, 1:3) = [0 0 1];
    P_add(n + i) = 0;
end

if ( strcmp(method, 'simple') )
    
    for i = 2:n-1
        beta_i = A(i - 1, i);
        alpha_i = A(i - 1, i - 1);
        
        d(i-1) = - alpha_i / beta_i;
        for j = 1 : p
            gamma_i_j = A(i - 1, n + j);
            e(i-1, j) = - gamma_i_j / beta_i;
        end
        f(i-1) = b(i - 1) / beta_i;
        
        D = 1;
        E = zeros(1, p);
        F = 0;
        for j = 1:i - 1
            D = D * d(j);
            
            tmp = 1;
            for k = j + 1 : i - 1
                tmp = tmp * d(k);
            end
            
            for k = 1 : p
                E(k) = E(k) + e(j, k) * tmp;
            end
            
            F = F + f(j) * tmp;
        end
        
        P(i, 1 : 2 + p) = [ D 0 E];
        P_add(i) = F;
    end
    
elseif (strcmp(method, 'fast'))
    
    D = 1;
    E = zeros(1, p);
    F = 0;
    
    for i = 2:n-1
        beta_i = A(i - 1, i);
        alpha_i = A(i - 1, i - 1);
        
        d(i-1) = - alpha_i / beta_i;
        for j = 1 : p
            gamma_i_j = A(i - 1, n + j);
            e(i-1, j) = - gamma_i_j / beta_i;
        end
        f(i-1) = b(i - 1) / beta_i;
        
        D = D * d(i-1); 
        F = f(i - 1) + d(i-1) * F;
        for j = 1:p
             E(j) = e(i - 1) + d(i-1) * E(j);
        end
        
        P(i, 1 : 2 + p) = [ D 0 E];
        P_add(i) = F;
    end
    
else
    throw (MException ('ArgumentCheck:IllegalArgument', 'Unsupported method!'));
end

A_ = A * P;
b_ = b - A * P_add;

A_ = A_(n - 2 : n-1, 1 : 2 + p);
b_ = b_(n - 2 : n-1);

% 
% A_ = zeros(2+)
% for i = 1:p+2
%         if ( i == 1 )
%             r
%         else
%             row = i + n-1;
%         end
% end


end

