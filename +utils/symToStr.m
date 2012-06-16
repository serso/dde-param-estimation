function result = symToStr( s )
%SYMTOSTR Summary of this function goes here
%   Detailed explanation goes here
result = char(s);

result = strrep(result, 'matrix', '');
result = substring(result, 2, length(result) - 2);
result = strrep(result, '[', '');
result = strrep(result, '],', ';');
result = strrep(result, ']]', '');

end

