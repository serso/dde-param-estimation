function reorderMatrix_TEST( M, b)
%%SOLVEMATRIX_TEST Summary of this function goes here
%   Detailed explanation goes here

methods = {'symrcm', 'amd', 'colamd', 'colperm', 'dmperm', 'symamd'};

for i = 1 : length(methods)
    method = methods(i);
    
    if ( strcmp(method, 'symrcm') )
        p = symrcm(M);
    elseif ( strcmp(method, 'amd') )
        p = amd(M);
    elseif ( strcmp(method, 'colamd') )
        p = colamd(M);
    elseif ( strcmp(method, 'colperm') )
        p = colperm(M);
    elseif ( strcmp(method, 'dmperm') )
        p = dmperm(M);
    elseif ( strcmp(method, 'symamd') )
        p = symamd(M);
    else
        throw (MException ('IllegalArgument:NotSupportedMethod', 'Method is not supported'));
    end
    
    Mp = M(p,p);
    
    subplot( ceil( length(methods) / 3 ), 3, i);
    
    spy(M, 'b');
    hold on;
    spy(Mp, 'r');
    
    title(method);
    legend('Before', 'After');
    
    bp = b(p);
    xp = Mp \ bp;
    [~, rp] = sort(p);
    x = xp(rp);
    
    solError = norm(M*x-b, inf);
    if  solError > 10^-10
        display(solError);
        throw (MException ('AssertionError:ConditionFailed', 'Not a solution!'));
    end
    
end

%%

end
