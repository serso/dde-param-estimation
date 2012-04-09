clc;

r = 10;
s = 5;

B = 2 * sparse(eye(r));
V = sprand (r, s, 0.2);
C_ = sprand (s, s, 0.2);

A = [B V; V' C_];

A = A'*A;

B = A(1:r, 1:r);
V = A(1:r, r:r+s);
C_ = A(r:r+s, r:r+s);

if ( norm(V'-A(r:r+s, 1:r), inf) > 0 ) 
    throw (MException ('AssertionError:Illegal', 'V != transpose(V) !'));
end

% spy(A);

[R, p] = chol(A);
if ( norm(p, inf) > 0 ) 
    throw (MException ('AssertionError:NotPositiveDefinite', 'Matrix is not positive definite!'));
end

Lb = chol(B, 'lower');
W = Lb \ V;
C = C_ - W'*W;
Lc = chol(C, 'lower');
L = [Lb zeros(r, s); W' Lc];