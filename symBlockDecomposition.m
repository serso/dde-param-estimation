function [ L, D ] = symBlockDecomposition( A, n )
%SYMBLOCKDECOMPOSITION decomposes block matrix of special form into product
%of 3 matrices A = L * D * L'. Matrix A must be of the next form:
% A = [ B V; V' 0]
% where 
% V' =
%     a1     b1     0      0      ...     0      g1_1       g2_1        gp_1
%     0      a2     b2     0      ...     0      g1_2       g2_2        gp_2
%     0      0      a3     b3     ...     0      g1_3       g2_3        gp_3
%     t1     0      0      a4     ...     0      g1_4       g2_4        gp_4
%     0      t2     0      0      ...     0      g1_5       g2_5        gp_5
%     ...    ...    ...    ...    ...     ...    ...        ...         ...
%     0      0      t(n-i) 0      0       b(n-1) g1_(n-1)   g2_(n_1)    gp_(n-1)
%
% B = 
%     b1     0      0      ...    0
%     0      b2     0      ...    0
%     0      0      b3     ...    0
%     ...    ...    ...    ...    0
%     0      0      0      ...    bn
%
%     bi >= 0
%
% Result matrices will be:
% D = 
%     d1     0      0      ...    0
%     0      d2     0      ...    0
%     0      0      d3     ...    0
%     ...    ...    ...    ...    0
%     0      0      0      ...    dn
%
% L = [ Lb 0; W' Lc ]

[rows, cols] = size(A);

% elements = zeros (n, 3);
% for i = 1:n
%     if ( abs(A(i, i)) > 0 )
%         [elements, ~] = addSparseElement(elements, i, i, i, 1);
%     else
%         [elements, ~] = addSparseElement(elements, i, i, i, 0);
%     end
% end
% Lb = createSparseMatrix(elements);
Lb = speye(n);
D1 = 2 * Lb;

%[Lb, D1, P] = ldl(A(1:n,1:n));
%Lb = P'\Lb;

V = A(1:n,n+1:cols);
W = V / 2;

method = 'asym';

if ( strcmp(method, 'sym') )
    tic;
    W0 = -W' * D1 * W;
    toc
    tic;
    [Lc, D2, P] = ldl(W0);
    toc;
elseif (strcmp(method, 'asym'))
    tic;
    W_ = Lb \ W;
    W0 = - V' * W_;
    toc
    tic;
    [Lc, D2, P] = ldl(W0);
    toc
else   
    throw (MException ('IllegalArgument:UnsupportedMathod', 'Unsupported method!'));
end

Lc = P' \  Lc;

O = zeros(n, cols-n);
L = [ Lb O; W' Lc ];
D = [D1 O; O' D2];

end

