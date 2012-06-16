function m = filterSmallElements( m, eps )
%FILTERSMALLELEMENTS removes (sets to 0) all elements in matrix which are less than eps in
%absolute value

[rows, cols] = size(m);
for row = 1:rows
    for col = 1:cols
        if ( abs(m(row, col)) <= eps )
            m(row, col) = 0;
        end
    end
end


end

