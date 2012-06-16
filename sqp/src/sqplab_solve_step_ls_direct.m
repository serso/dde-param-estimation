function [ dlmp, p, rp ] = sqplab_solve_step_ls_direct( M, me, info, options )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

     A = [M, info.ae'; info.ae, sparse(me, me)];
     b = [-info.g; -info.ce];
     
     p = [];
     rp = [];
     
     if ( strcmp(options.stepMethod, 'reduce') )
         
         [P, P_add, ae_, ce_] = p_factors(info.ae, -info.ce, 1);
         
         M_ = transform_m(M, 1);
         g_ = transform_m(info.g, 1);
         
         [rows, ~] = size(ae_);
         Z = sparse(rows, rows);
         
         A_ = [M_, ae_'; ae_, Z];
         %           spy(A_);
         
         %           dlmp_ = A_ \ [-g_; ce_];
         dlmp_ = A_ \ [-g_; ce_];
         pk_ = dlmp_(1:3);
         
         dlmp = P * pk_ + P_add;
         
     elseif ( strcmp(options.stepMethod, 'ldls') )
         if ( isfield(options, 'ldlsThreshold') )
             [L,D,P,S] = ldl(A, options.ldlsThreshold);
         else
             [L,D,P,S] = ldl(A);
         end
         dlmp = (S*P) * (L'\(D\(L\( (P'*S) * b))));
     elseif (strcmp(options.stepMethod, 'block-decomposition'))
         [n, ~] = size(M);
         [L, D] = blockDecomposition(A, n);
         dlmp = L'\(D\(L\b));
     elseif (strcmp(options.stepMethod, 'ldl'))
         %tic;
         [L,D,P] = ldl(A);
         dlmp= P * (L'\(D\(L\(P' * b))));
         %toc
     elseif (strcmp(options.stepMethod, 'lu'))
         [L, U, P, Q] = lu(A);
         y = L \ (P * b);
         dlmp = Q * (U \ y);
     elseif ( strcmp(options.stepMethod, 'symrcm') || ...
             strcmp(options.stepMethod, 'amd') || ...
             strcmp(options.stepMethod, 'colamd') || ...
             strcmp(options.stepMethod, 'colperm') || ...
             strcmp(options.stepMethod, 'dmperm') || ...
             strcmp(options.stepMethod, 'symamd'))
         if ( ~isfield(options, 'p') )
             options.p= [];
         end
         
         if (  isempty(options.p) )
             
             if ( strcmp(options.stepMethod, 'symrcm') )
                 options.p = symrcm(A);
             elseif ( strcmp(options.stepMethod, 'amd') )
                 options.p = amd(A);
             elseif ( strcmp(options.stepMethod, 'colamd') )
                 options.p = colamd(A);
             elseif ( strcmp(options.stepMethod, 'colperm') )
                 options.p = colperm(A);
             elseif ( strcmp(options.stepMethod, 'dmperm') )
                 options.p = dmperm(A);
             elseif ( strcmp(options.stepMethod, 'symamd') )
                 options.p = symamd(A);
             else
                 throw (MException ('IllegalArgument:NotSupportedMethod', 'Method is not supported'));
             end
             [~, rp] = sort(options.p);
             options.rp = rp;
         end
         
         p = options.p;
         rp = options.rp;
         
         if ( strcmp(options.stepMethod, 'symrcm') )
             
             [L, U, P] = lu(A(p,p));
             y = L \ (P * b(p));
             x0p = U \ y;
             dlmp = x0p(rp);
             
         elseif ( strcmp(options.stepMethod, 'amd') )
             % for Cholecky only
             throw (MException ('IllegalArgument:NotSupportedMethod', 'Method is not supported'));
         elseif ( strcmp(options.stepMethod, 'colamd') )
             
             [L, U, P] = lu(A(:,p));
             y = L \ (P * b);
             x0p = U \ y;
             dlmp = x0p(rp);
             
         elseif ( strcmp(options.stepMethod, 'colperm') )
             
             [L, U, P] = lu(A(:,p));
             y = L \ (P * b);
             x0p = U \ y;
             dlmp = x0p(rp);
             
         elseif ( strcmp(options.stepMethod, 'dmperm') )
             % nothing useful in documentation found
             throw (MException ('IllegalArgument:NotSupportedMethod', 'Method is not supported'));
         elseif ( strcmp(options.stepMethod, 'symamd') )

             [L, U, P] = lu(A(p,p));
             y = L \ (P * b(p));
             x0p = U \ y;
             dlmp = x0p(rp);
             
         else
             throw (MException ('IllegalArgument:NotSupportedMethod', 'Method is not supported'));
         end
     elseif ( strcmp(options.stepMethod, 'mldivide') )
         dlmp =  A \ b;
     else
         display(options.stepMethod);
         throw (MException ('IllegalArgument:NotSupportedMethod', 'Method is not supported'));
     end
     
end

