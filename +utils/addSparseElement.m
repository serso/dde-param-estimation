%%
% Method adds value to elements matrix (where elements matrix is input matrix for createSparseMatrix() method)
%
% @param elements matrix where value will be added (must be Nx3)
% @param index sequential number - position in elements matrix where
% current entry will be stored (method affects this value)
% @param i - number of row where value will be stored
% @param j - number of column where value will be stored
% @param value - value to be stored in sparse matrix in the (i, j) position
function [elements, indexResult] = addSparseElement(elements, index, i, j, value)
    elements(index, :) = [i j value];
    indexResult = index + 1;
end
