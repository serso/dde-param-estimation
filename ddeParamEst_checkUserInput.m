function [ f, fg, fh, fSym, fgSym, fhSym, xHistory, xHistorySym ] = ddeParamEst_checkUserInput( f, np, delays, userConstants, xHistory )
%DDEPARAMEST_CHECKUSERINPUT Function checks user input data and convert it
%to MATLAB objects if no error occurred;

%%
% first check input data

if ( np < 1 )
    throw (MException ('ArgumentCheck:IllegalArgument', 'Number of estimated parameters "np" must be more than 0!'));
end

for i=1:np
    p_i = strcat('p_', num2str(i));
    if (isempty(strfind(f, p_i)))
        warning(p_i);
        throw (MException ('ArgumentCheck:IllegalArgument', 'Parameter p_i is not present in input function string!'));
    end
end

%%
% then create matlab objects
t = sym('t');

clearvars fUknowns;
fVars = cell(1, 3);

x = sym('x');
fUnknowns(1) = x;

realDelays = 0;
for i = 1:delays
    x_i = strcat('x_', num2str(i));
    if (~isempty(strfind(f, x_i)))
        fUnknowns(2 + realDelays) = sym(x_i);
        realDelays = realDelays + 1;
    end
end

fVars{1} = transpose(fUnknowns);
fVars{2} = t;

for i = 1:np
    fUnknowns(i + 1 + realDelays) = sym(strcat('p_', num2str(i)));
end

fVars{3} = transpose(fUnknowns(2 + realDelays:end));

fSym = evalin(symengine, f);

for i =  1 : length(userConstants)
    userConstant = userConstants(i);
    if ( ~isempty(userConstant{1}{1}) && ~isempty(userConstant{1}{2}) )
        fSym = subs(fSym, userConstant{1}{1}, userConstant{1}{2}, 0);
    end
end

fgSym = jacobian(fSym, fUnknowns);
fhSym = jacobian(fgSym, fUnknowns);

f = matlabFunction(fSym, 'vars', fVars);
fg = matlabFunction(fgSym, 'vars', fVars);
fh = matlabFunction(fhSym, 'vars', fVars);

%% xHistory

if ( nargin > 4 )
    xHistorySym = evalin(symengine, xHistory);
    
    for i =  1 : length(userConstants)
        userConstant = userConstants(i);
        if ( ~isempty(userConstant{1}{1}) && ~isempty(userConstant{1}{2}) )
            xHistorySym = subs(xHistorySym, userConstant{1}{1}, userConstant{1}{2}, 0);
        end
    end
    
    t = sym('t');
    
    xHistory = matlabFunction(xHistorySym, 'vars', {t});
else
    xHistorySym = [];
    xHistory = [];
end

end

