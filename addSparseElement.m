function [elements, indexResult] = addSparseElement(elements, index, i, j, value)
    elements(index, :) = [i j value];
    indexResult = index + 1;
end
