function [ sparseElements, i ] = symSparseElements( sparseElements, i )
%SYMSPARSEELEMENTS Summary of this function goes here
%   Detailed explanation goes here
% do array symmetric
for sparseElement = sparseElements'
    % i < j
    if ( sparseElement(2) < sparseElement(1) )
        [sparseElements, i] = utils.addSparseElement(sparseElements, i, sparseElement(2), sparseElement(1), sparseElement(3));
    end
end

end

