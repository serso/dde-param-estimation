function mResult = createSparseMatrix(elements)
    mResult = sparse(elements( : , 1), elements( : , 2), elements( : , 3));
end