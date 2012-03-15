function [ a_ ] = transform_m( a, p )
%TRANSFORM_M Summary of this function goes here
%   Detailed explanation goes here

[rows, cols] = size(a);
if ( rows == 1 )
    a_ = zeros(1, 2 + p);    
    
    a_ (1, 1) = a (1, 1);
    for i = cols - p: cols
        a_ (1, i - cols + 2 + p) = a (1, i);
    end
    
elseif ( cols == 1 )
    a_ = zeros(2 + p, 1);    
    a_ (1) = a (1);
    for i = rows - p: rows
        a_ (i - rows + 2 + p) = a (i);
    end
else
    a_ = zeros(2 + p, 2 + p); 
    
    a_(1, 1) = a(1, 1);
    for i = 2: 2 + p
        a_(i, 1) = a(i + rows - 2 - p, 1);
        a_(1, i) = a(1, i + cols - 2 - p);
    end
    
    for i = rows - p: rows
        for j = cols - p: cols
            a_ (i - rows + 2 + p, j - cols + 2 + p) = a (i, j);
        end
    end

end;

end

