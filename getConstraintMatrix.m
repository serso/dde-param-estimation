function [ce] = getConstraintMatrix ( t, N, p, delays, delayF, method)

ceElements = zeros(2 * (N - 1) + p * (N - 1), 3);
ceElementI = 1;

if ( strcmp(method, 'euler') )
    
    for i = 1:1:(N - 1)
        delta_t_i = delta(t, i);
        
        [ceElements, ceElementI] = addSparseElement(ceElements, ceElementI, i, i, delta_t_i - 1);
        [ceElements, ceElementI] = addSparseElement(ceElements, ceElementI, i, i+1, 1);
        
        for j = (N+1):1:(N+p)
            [ceElements, ceElementI] = addSparseElement(ceElements, ceElementI, i, j, -delta_t_i);
        end
    end
    
elseif (strcmp(method, 'backward_euler'))
    
    for i = 1:1:(N - 1)
        delta_t_i = delta(t, i);
        
        [ceElements, ceElementI] = addSparseElement(ceElements, ceElementI, i, i, -1);
        [ceElements, ceElementI] = addSparseElement(ceElements, ceElementI, i, i+1, 1 + delta_t_i );
        
        for j = (N+1):1:(N+p)
            [ceElements, ceElementI] = addSparseElement(ceElements, ceElementI, i, j, -delta_t_i );
        end
    end
    
elseif (strcmp(method, 'box'))
    
    for i = 1:1:(N - 1)
        delta_t_i = delta(t, i);
        
        [ceElements, ceElementI] = addSparseElement(ceElements, ceElementI, i, i, 1 - delta_t_i / 2);
        [ceElements, ceElementI] = addSparseElement(ceElements, ceElementI, i, i+1, 1 + delta_t_i / 2 );
        
        for j = (N+1):1:(N+p)
            [ceElements, ceElementI] = addSparseElement(ceElements, ceElementI, i, j, -delta_t_i);
        end
    end
    
elseif (strcmp(method, 'rk4'))
    
    
else
    throw (MException ('ArgumentCheck:IllegalArgument', 'Unsupported method!'));
end

ce = createSparseMatrix(ceElements);

% if ( ~isempty(ce) )
%     spy(ce);
% end

end