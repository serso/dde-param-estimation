function [ s ] = setDefaultOption( s, fieldName, fieldValue )
%DDEPARAMEST_SETDEFAULTOPTION Sets default field value if no field value is
%set
%s = DDEPARAMEST_SETDEFAULTOPTION(s, fieldName, fieldValue) checks is field
% with name 'fieldName' is present in the structure 's' and sets 'fieldValue' if
% no such field exists
if ( ~isfield(s, fieldName) )
    s.(fieldName) = fieldValue;
end

end

