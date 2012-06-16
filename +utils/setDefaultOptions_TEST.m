function setDefaultOptions_TEST
%%
clc;
s = [];
s.('test4') = 12;
s = utils.setDefaultOptions(s, {{'test', true}, {'test2', false}, {'test3', 45}, {'test4', 42}});
if ( ~isfield(s, 'test') )
    throw (MException('AssertionError:ConditionFaild', 'test is not a field of s')) ;
end

if ( ~s.('test') )
     throw (MException('AssertionError:ConditionFaild', 'Incorrect value for test field')) ;
end

if ( ~isfield(s, 'test2') )
    throw (MException('AssertionError:ConditionFaild', 'test2 is not a field of s')) ;
end

if ( s.('test2') )
     throw (MException('AssertionError:ConditionFaild', 'Incorrect value for test2 field')) ;
end

if ( ~isfield(s, 'test3') )
    throw (MException('AssertionError:ConditionFaild', 'test3 is not a field of s')) ;
end

if ( s.('test3') ~= 45 )
     throw (MException('AssertionError:ConditionFaild', 'Incorrect value for test3 field')) ;
end

if ( s.('test4') ~= 12 )
     throw (MException('AssertionError:ConditionFaild', 'Incorrect value for test4 field')) ;
end

%%
clc;
s = [];
s.('test4') = 12;
s = utils.setDefaultOptions(s, ...
    {
        {'test', true}
        {'test2', false}
        {'test3', 45}
        {'test4', 42}
    });
if ( ~isfield(s, 'test') )
    throw (MException('AssertionError:ConditionFaild', 'test is not a field of s')) ;
end

if ( ~s.('test') )
     throw (MException('AssertionError:ConditionFaild', 'Incorrect value for test field')) ;
end

if ( ~isfield(s, 'test2') )
    throw (MException('AssertionError:ConditionFaild', 'test2 is not a field of s')) ;
end

if ( s.('test2') )
     throw (MException('AssertionError:ConditionFaild', 'Incorrect value for test2 field')) ;
end

if ( ~isfield(s, 'test3') )
    throw (MException('AssertionError:ConditionFaild', 'test3 is not a field of s')) ;
end

if ( s.('test3') ~= 45 )
     throw (MException('AssertionError:ConditionFaild', 'Incorrect value for test3 field')) ;
end

if ( s.('test4') ~= 12 )
     throw (MException('AssertionError:ConditionFaild', 'Incorrect value for test4 field')) ;
end

%%
end

