
function setDefaultOption_TEST
%%

s.('test2') = true;
s = utils.setDefaultOption(s, 'test', false);
s = utils.setDefaultOption(s, 'test2', false);
if ( ~isfield(s, 'test') )
    throw (MException('AssertionError:ConditionFaild', 'test is not a field of s')) ;
end

if ( ~s.('test2') )
    throw (MException('AssertionError:ConditionFaild', 'test2 is not true of s')) ;
end



%%
end

