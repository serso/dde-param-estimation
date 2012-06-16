function [ s ] = setDefaultOptions( s, defaultFields )
%DDEPARAMEST_SETDEFAULTOPTIONS Sets default field values for all fields
%which are not set
%s = DDEPARAMEST_SETDEFAULTOPTIONS(s, defaultFields) iterates over
%'defaultFields' array and for each 'defaultField' checks if field present
%and set default value. Each 'defaultField' must be a cell array of 2
%elements
for i = 1:length(defaultFields)
    defaultField = defaultFields(i);
    s = utils.setDefaultOption(s, defaultField{1}{1}, defaultField{1}{2});
end

end

