%%
% Method creates sparse matrix from specified 'elements' matrix.
%
% @param elements - Nx3 matrix, where
% first column contains number of result matrix row (i)
% second column contains number of result matrix column (j)
% third column contains value in i-th row and j-th column of result matrix.
%
% e.g. elements =   [1, 2, 10.0;
%                    2, 3, 15.0]
% yields matrix 3-by-3 where all elements are zeroes except (1,2)=10.0 and
% (2,3)=15.0
function mResult = createSparseMatrix(elements)
    mResult = sparse(elements( : , 1), elements( : , 2), elements( : , 3));
end