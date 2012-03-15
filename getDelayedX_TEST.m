function getDelayedX_TEST()
%%GETDELAYEDX_TEST Summary of this function goes here
%   Detailed explanation goes here


% Test 01

x = 1:10;
t = 1:10;
delays = [3 10];
delayF = @(t) t;

delayedX = getDelayedX(x, t, 5, delays, delayF);

if (norm(delayedX - [2; -5], inf) > 0)
     throw (MException ('AssertionError:ConditionFailed', 'Expected result is different than actual!')); 
end

% Test 02
x = 1:0.1:10;
t = 1:0.1:10;
delays = [3 10];
delayF = @(t) t;

delayedX = getDelayedX(x, t, 41, delays, delayF);
if (norm(delayedX - [2; -5], inf) > 0)
     throw (MException ('AssertionError:ConditionFailed', 'Expected result is different than actual!')); 
end

% Test 03
x = -10:0.1:-1;
t = -10:0.1:-1;
delays = [3 10];
delayF = @(t) t;

delayedX = getDelayedX(x, t, 41, delays, delayF);
if (norm(delayedX - [-9; -16], inf) > 0)
     throw (MException ('AssertionError:ConditionFailed', 'Expected result is different than actual!')); 
end

%%
end

