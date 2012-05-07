function [xResult, lambdaResult, iterations, funCount, timeResult, sqpInfo] = sqp ( N, cn, f, lh, ce, cej, x0, options, debug, lb, ub)

%
% [xResult, lambdaResult, iterations, funCount, xStepsResult, timeResult] = sqp (N, cn, f, l_grad, lh, ce, cej, x0, options,
% debug, lb, ub);
%

timerId = tic();

if ( isempty(N) )
    throw (MException ('ArgumentCheck:NotSet', 'N is not set.'));
elseif ( N <= 0 )
    throw (MException ('ArgumentCheck:OutOfRange', 'N is less than 1.'));
end

if ( isempty(lh) )
    throw (MException ('ArgumentCheck:NotSet', 'Hessian must be set.'));
end

if ( isempty(x0) )
    x0 = zeros(N, 1);
elseif ( length(x0) ~= N )
    throw (MException ('ArgumentCheck:illegalArgument', 'x0 length is not equal to N.'));
end

iterations = 0;
funCount = 0;

    function [outdic, out1, out2, out3, out4, out5, out6, out7  ] = dde_simul ( indic, x, lm )
        
        out1 = [];
        out2 = [];
        out3 = [];
        out4 = [];
        out5 = [];
        out6 = [];
        out7 = [];
        
        % timerId = tic();
        
        % display(indic);
        
        if indic <= 1
            outdic = 0;
            iterations = iterations + 1;
        elseif indic == 2
            % f,ci,ce,cs,g,ai,cej
            outdic = 0;
            out1 = f(x);
            funCount = funCount + 1;
            out3 = ce(x);
        elseif indic == 3
            % f,ci,ce,cs,g,ai,cej
            outdic = 0;
            funCount = funCount + 1;
            [~, out5] = f(x);
            out7 = cej(x);
        elseif indic == 4
            % f,ci,ce,cs,g,ai,cej
            outdic = 0;
            
            funCount = funCount + 1;
            [out1, out5] = f(x);
            out3 = ce(x);
            out7 = cej(x);
        elseif indic == 5
            outdic = 0;
            out1 = lh(x, lm);
        else
            outdic = -2;
        end
        
        %time = toc(timerId);
        %display(sprintf('indic %i %0.5f s', indic, time));
    end

[xResult,lambdaResult,sqpInfo] = sqplab(@dde_simul, x0, zeros(N + cn, 1), lb, ub, options);

timeResult = toc(timerId);

if ( debug )
    display(sprintf('Time: %0.5f s', timeResult));
end

end