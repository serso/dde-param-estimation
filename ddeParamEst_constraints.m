function [ci, ce, cij, cej, cih, ceh] = ddeParamEst_constraints ( ...
    f, fg, fh,      ...
    t, x,           ...
    n, np,           ...
    delays, xHistory,      ...
    method, options, lm)
%DDEPARAMEST_COnSTRAInTS calculates constraints in point x for the DDE parameters
%estimation problem: (*) dx/dt = f (t, x(t), x(t-τ_1), ..., x(t - τ_n))
%
%   [ci, ce, cij, cej] = ddeParamEst_constraints(f, fg, fh, t, x, n, np, delays, xHistory, method, lm) calculates next values:
%   ci  vector of inequalities calculated in the point (t, x)
%   ce  vector of equalitites calculated in the point (t, x)
%   cij jacobian matrix of inequalities in the point (t, x)
%   cej jacobian matrix of equalities in the point (t, x)
%   cih hessian matrix of inequalities in the point (t, x)
%   ceh hessian matrix of equalities in the point (t, x)
%
%
% @param f - handler to the function f of (*) (might be [])
% @param fg - handler to the gradient of f (might be [])
% @param fh - handler to the hessian of f (might be [])
% @param t - time grid
% @param x - x grid
% @param n - number of elements in integration grid
% @param np - number of estimated parameters
% @param delays - vector of delays tau (if delay is estimated - use nUMBER OF ELEMEnT in parameter estimation vectors )
% @param xHistory - function for history of x in [tMin - max(τ_i), tMin]
% @param method - derivative approximation method
%

%% INIT
% prepare initial parameters

% extract parameters
p = x ( n + 1 : n + np );

% transform delays (substitute estimates delays with values from p)
origDelays = delays;
delays = getDelays(delays, p);

ndelays = length(delays);

% if no f handler is provided => no need for calculations of ce => set empty
if (~isempty(f))
    % allocate ce array
    ce = zeros (n - 1, 1);
else
    ce = [];
end

% prevailing index for cej
ceji = 1;
% array of values of jacobian
ceja = [];
cejaLength = 0;

% if no fg handler is provided => no need for calculations of cej => set
% empty
if ( ~isempty(fg) )
    % count total value of non zero element to effectively allocate memory
    delayIndeces = getDelayIndeces(t, delays, p, options);
    
    % number of elements in the "diagonal"
    ndiag = n-1;
    
    totalDelayElements = 0;
    for delayIndex = delayIndeces
        
        delayElements = 0;
        
        if (delayIndex > 1)
            % NOTE: diagonal is already counted in the code below
            
            % each delay has a sub diagonal => number of not zero elements
            % can be retrieved by the next formula
            if ( strcmp(method, 'euler') )
                delayElements = (ndiag - (delayIndex - 1));
            elseif (strcmp(method, 'backward-euler'))
%                 if ( delayIndex == 2 )
%                     display(delayIndex)
%                     throw (MException ('ArgumentCheck:IllegalArgument', 'Check!'));
%                 end
                delayElements = (ndiag - (delayIndex - 1) + 1);
            elseif (strcmp(method, 'box'))
                delayElements = (n - delayIndex - 2);
            else
                throw (MException ('ArgumentCheck:IllegalArgument', 'Unsupported method!'));
            end
        end
        
        if ( delayElements > 0 )
            totalDelayElements = totalDelayElements + delayElements;
        end
        
    end
    
    % allocating memory
    ceja = zeros(2 * ndiag + np * (n - 1) + totalDelayElements, 3);
    
    cejaLength = length(ceja);
else
    cej = [];
end

% prevailing index for ceh
cehi = 1;

% array of values of hessian
ceha = [];

cehaLength = 0;

if ( ~isempty(fh) )
    % count total value of non zero element to effectively allocate memory
    delayIndeces = getDelayIndeces(t, delays, p, options);
    
    % matrix of % d²(l'c)/d(pi)d(pj)
    dpi_dpj = zeros(np);
    
    % number of elements in the "diagonal"
    ndiag = n;
    
    totalDelayElements = 0;
    
    for delayIndex = delayIndeces
        
        delayElements = 0;
        
        if (delayIndex > 1)
            % NOTE: diagonal is already counted in the code below
            
            % each delay has a sub diagonal => number of not zero elements
            % can be retrieved by the next formula
            if ( strcmp(method, 'euler') )
                delayElements = ndiag - (delayIndex - 1);
            elseif (strcmp(method, 'backward-euler'))
                delayElements = (ndiag - (delayIndex - 1)) - 1;
            elseif (strcmp(method, 'box'))
                throw (MException ('ArgumentCheck:IllegalArgument', 'Unsupported method!'));
            else
                throw (MException ('ArgumentCheck:IllegalArgument', 'Unsupported method!'));
            end
        end
        
        if ( delayElements > 0 ) 
            totalDelayElements = totalDelayElements + delayElements;
        end
    end
       
    % don't forget symmetric elements
    % and mixed derivatives
    totalDelayElements = 2 * totalDelayElements;
    
    % allocating memory
    ceha = zeros((n - 1) + 2 * np * (n - 1) + np * np + totalDelayElements, 3);
    
    cehaLength = length(ceha);
    
    if ( isempty(lm) )
        lm = ones(n-1, 1);
    end
else
    ceh = [];
end

%% CALCULATE
% create output matrices

if ( strcmp(method, 'euler') )
    
    %% EULER
    for i = 1:(n - 1)
        
        delayedX = getDelayedX(x, t, i, delays, xHistory, options);
        
        if ( ~isempty(f) )
            ce(i) = - x(i) + x(i + 1) - delta(t, i) * f(delayedX, t(i), p);
        end
        
        if ( ~isempty(fg) )
            
            [xd, pd] = getDerivatives(fg(delayedX, t(i), p), ndelays, np);
            
            [ceja, ceji] = utils.addSparseElement(ceja, ceji, i, i, (- 1 - delta(t, i) * xd(1)));
            [ceja, ceji] = utils.addSparseElement(ceja, ceji, i, i+1, 1);
            
            for pi = 1:np
                [ceja, ceji] = utils.addSparseElement(ceja, ceji, i, n + pi, (- delta(t, i) * pd(pi)));
            end
            
            xdIndex = 1;
            for delayIndex = delayIndeces
                if delayIndex > 1
                    
                    % delay index relative to current position. Note: +1 as delay index for tau = 0
                    % is 1
                    di = i - (delayIndex - 1);
                    if ( di > 0 )
                        [ceja, ceji] = utils.addSparseElement(ceja, ceji, i, di, - delta(t, i) * xd(xdIndex));
                    end
                end
                xdIndex = xdIndex + 1;
            end
        end
        
        if (~isempty(fh))            
            % fhi - function hessian at the point (t_i, x_i)
            fhi = deal(fh(getDelayedX(x, t, i, delays, xHistory, options), t(i), p));
            
            % d²(l'c)/d(xi²)
            dxi_dxi= - lm(i) * delta(t, i) * fhi(1, 1);
            
            % d²(l'c)/d(xi)d(pi)
            dxi_dpj = zeros(np, 1);
            for pi = 1:np
                dxi_dpj(pi) = - lm(i) * delta(t, i) * fhi(1, ndelays + pi);
            end
            
            fhiIndex = 1;
            for delayIndex = delayIndeces
                if delayIndex > 1
                    
                    j = i + (delayIndex - 2);
                    if ( j < n )
                        % fhdi - hessian as the point (t_delayIndex, x_delayIndex)
                        fhdi = deal(fh(getDelayedX(x, t, j, delays, xHistory, options), t(j), p));
                        
                        % diagonal elements
                        dxi_dxi = dxi_dxi - lm(j) * delta(t, j) * fhdi(1, 1);
                        
                        % delay elements
                        %[ceha, cehi] = utils.addSparseElement(ceha, cehi, j, i, - lm(i) * delta(t, i) * fhi(1, fhiIndex));
                        
                        % parameter elements
                        for pi = 1:np
                            dxi_dpj(pi) = dxi_dpj(pi) - lm(j) * delta(t, j) * fhdi(fhiIndex, ndelays + pi);
                        end
                    end
                    
                    j = i - (delayIndex - 2);
                    if ( j > 0 )
                        [ceha, cehi] = utils.addSparseElement(ceha, cehi, i, j, - lm(i) * delta(t, i) * fhi(1, fhiIndex));
                    end
                end
                fhiIndex = fhiIndex + 1;
            end
            
            [ceha, cehi] = utils.addSparseElement(ceha, cehi, i, i, dxi_dxi);
            
            for pi = 1:np
                [ceha, cehi] = utils.addSparseElement(ceha, cehi, n + pi, i, dxi_dpj(pi));
            end
            
            for pi = 1:np
                for pj = 1:np
                    dpi_dpj(pi, pj) = dpi_dpj(pi, pj) - lm(i) * delta(t, i) * fhi(ndelays + pi, ndelays + pj);
                end
            end
        end
        
    end
    
    if (~isempty(fh))
        for pi = 1:np
            for pj = 1:pi
                [ceha, cehi] = utils.addSparseElement(ceha, cehi, n + pi, n + pj, dpi_dpj(pi, pj));
            end
        end
    end
    
    if ( ~isempty(fg) )
        cej = utils.createSparseMatrix(ceja);
    end
    
    if ( ~isempty(fh) )
        ceh = utils.createSparseMatrix(utils.symSparseElements(ceha, cehi));
    end
    
elseif (strcmp(method, 'backward-euler'))
    
    %% BACKWARD-EULER
    for i = 1:(n - 1)
        
        delayedX = getDelayedX(x, t, i + 1, delays, xHistory, options);
        
        if ( ~isempty(f) )
            ce(i) = - x(i) + x(i + 1) - delta(t, i) * f ( delayedX, t ( i + 1), p );
        end;
        
        if ( ~isempty(fg) )
            
            [xd, pd] = getDerivatives(fg(delayedX, t(i+1), p), ndelays, np);
            
            [ceja, ceji] = utils.addSparseElement(ceja, ceji, i, i, -1);
            [ceja, ceji] = utils.addSparseElement(ceja, ceji, i, i+1, 1 - delta(t, i) * xd(1) );
            
            for pi = 1:np
                [ceja, ceji] = utils.addSparseElement(ceja, ceji, i, n + pi, -delta(t, i) * pd(pi));
            end
            
            xdIndex = 1;
            for delayIndex = delayIndeces
                if delayIndex > 1
                    
                    % delay index relative to current position. Note: +1 as delay index for tau = 0
                    % is 1
                    di = i - (delayIndex - 1) + 1;
                    if ( di > 0 )
                        [ceja, ceji] = utils.addSparseElement(ceja, ceji, i, di, - delta(t, i) * xd(xdIndex));
                    end
                end
                xdIndex = xdIndex + 1;
            end
        end
        
        if (~isempty(fh))
            % fhi - function hessian at the point (t_i, x_i)
            fhi = deal(fh(getDelayedX(x, t, i + 1, delays, xHistory, options), t(i + 1), p));
            
            % d²(l'c)/d(xi+1²)
            dxi1_dxi1 = - lm(i) * delta(t, i) * fhi(1, 1);
            
            % d²(l'c)/d(xi+1)d(pi)
            dxi1_dpj = zeros(np, 1);            
            for pi = 1:np
                dxi1_dpj(pi) = - lm(i) * delta(t, i) * fhi(1, ndelays + pi);
            end
            
            fhiIndex = 1;
            for delayIndex = delayIndeces
                if delayIndex > 1
                    
                    j = i + (delayIndex - 1);
                    if ( j < n )
                        % fhdi - hessian as the point (t_delayIndex, x_delayIndex)
                        fhdi = deal(fh(getDelayedX(x, t, j + 1, delays, xHistory, options), t(j + 1), p));
                        % diagonal elements
                        dxi1_dxi1 = dxi1_dxi1 - lm(j) * delta(t, j) * fhdi(1, 1);
                        
                        % delay elements
                        %[ceha, cehi] = utils.addSparseElement(ceha, cehi, j, i, - lm(i) * delta(t, i) * fhi(1, fhiIndex));
                        
                        % parameter elements
                        for pi = 1:np
                            dxi1_dpj(pi) = dxi1_dpj(pi) - lm(j) * delta(t, j) * fhdi(fhiIndex, ndelays + pi);
                        end
                    end
                    
                    j = i - (delayIndex - 1);
                    if ( j >= 0 )
                        [ceha, cehi] = utils.addSparseElement(ceha, cehi, i + 1, j + 1, - lm(i) * delta(t, i) * fhi(1, fhiIndex));
                    end
                end
                fhiIndex = fhiIndex + 1;
            end
            
            [ceha, cehi] = utils.addSparseElement(ceha, cehi, i+1, i+1, dxi1_dxi1);
            
            for pi = 1:np
                [ceha, cehi] = utils.addSparseElement(ceha, cehi, n + pi, i + 1, dxi1_dpj(pi));
            end
            
            for pi = 1:np
                for pj = 1:np
                    dpi_dpj(pi, pj) = dpi_dpj(pi, pj) - lm(i) * delta(t, i) * fhi(ndelays + pi, ndelays + pj);
                end
            end
        end
    end
    
    if (~isempty(fh))
        for pi = 1:np
            for pj = 1:pi
                [ceha, cehi] = utils.addSparseElement(ceha, cehi, n + pi, n + pj, dpi_dpj(pi, pj));
            end
        end
        
        dxi1_dpj = zeros(np, 1);
        dxi1_dxi1 = 0;
        
        % todo serso: actually zero step can be done inside main loop, indexes must be
        % changed.
        
        % do zero step
        i = 0;
        fhiIndex = 1;
        for delayIndex = delayIndeces
            if delayIndex > 1
                
                j = i + (delayIndex - 1);
                if ( j < n )
                    % fhdi - hessian as the point (t_delayIndex, x_delayIndex)
                    fhdi = deal(fh(getDelayedX(x, t, j + 1, delays, xHistory, options), t(j + 1), p));
                    % diagonal elements
                    dxi1_dxi1 = dxi1_dxi1 - lm(j) * delta(t, j) * fhdi(1, 1);
                    
                    % delay elements
                    %[ceha, cehi] = utils.addSparseElement(ceha, cehi, j, i, - lm(i) * delta(t, i) * fhi(1, fhiIndex));
                    
                    % parameter elements
                    for pi = 1:np
                        dxi1_dpj(pi) = dxi1_dpj(pi) - lm(j) * delta(t, j) * fhdi(fhiIndex, ndelays + pi);
                    end
                end
                
                j = i - (delayIndex - 1);
                if ( j >= 0 )
                    [ceha, cehi] = utils.addSparseElement(ceha, cehi, i + 1, j + 1, - lm(i) * delta(t, i) * fhi(1, fhiIndex));
                end
            end
            fhiIndex = fhiIndex + 1;
        end
        
        [ceha, cehi] = utils.addSparseElement(ceha, cehi, i+1, i+1, dxi1_dxi1);
        
        for pi = 1:np
            [ceha, cehi] = utils.addSparseElement(ceha, cehi, n + pi, i + 1, dxi1_dpj(pi));
        end
    end
    
    if ( ~isempty(fg) )
        cej = utils.createSparseMatrix(ceja);
    end
    
    if ( ~isempty(fh) )
        ceh = utils.createSparseMatrix(utils.symSparseElements(ceha, cehi));
    end
    
elseif (strcmp(method, 'box'))
    
    for i = 1:1:(n - 1)
        delayedX_i = getDelayedX(x, t, i, delays, xHistory, options);
        delayedX_i_1 = getDelayedX(x, t, i + 1, delays, xHistory, options);
        delayedX = (delayedX_i + delayedX_i_1) / 2;
        
        if ( i < n - 1)
            deltaT = (delta(t, i) + delta(t, i+1)) / 2;
        else
            deltaT = delta(t, i);
        end
        
        %deltaT = delta(t, i);
        %ti = t(i) + (t(i+1) - t(i)) / 2;
        ti = t(i) + deltaT;
        
        if ( ~isempty(f) )
            ce(i) = - x(i) + x(i + 1) - deltaT * f ( delayedX, ti, p );
        end
        
        if ( ~isempty(fg) )
            [xd, pd] = getDerivatives(fg(delayedX, ti, p), ndelays, np);
            
            [ceja, ceji] = utils.addSparseElement(ceja, ceji, i, i, ( - 1 - deltaT * xd(1) / 2) );
            [ceja, ceji] = utils.addSparseElement(ceja, ceji, i, i+1, (1 - deltaT * xd(1) / 2) );
            
            for j = (n+1):1:(n+np)
                [ceja, ceji] = utils.addSparseElement(ceja, ceji, i, j, (- deltaT * pd(j - n)) );
            end
            
            xdIndex = 1;
            for delayIndex = delayIndeces
                if delayIndex > 1
                    col = i - delayIndex;
                    if ( col > 0 )
                        [ceja, ceji] = utils.addSparseElement(ceja, ceji, i, col, - deltaT * xd(xdIndex) / 2);
                        [ceja, ceji] = utils.addSparseElement(ceja, ceji, i, col+1, - deltaT * xd(xdIndex) / 2);
                    end
                end
                xdIndex = xdIndex + 1;
            end
        end
        
    end
    
    if ( ~isempty(fg) )
        cej = utils.createSparseMatrix(ceja);
    end
    
elseif (strcmp(method, 'rk4'))
    
    cej = [];
    
    for i = 1:1:(n - 1)
        delayedX = getDelayedX(x, t, i, delays, xHistory, options);
        
        if ( ~isempty(f) )
            K1 = f(delayedX, t(i), p);
            K2 = f(delayedX + delta(t, i) * K1 / 2, t(i) + delta(t, i) / 2, p);
            K3 = f(delayedX + delta(t, i) * K2 / 2, t(i) + delta(t, i) / 2, p);
            K4 = f(delayedX + delta(t, i) * K3, t(i) + delta(t, i), p);
            
            ce(i) = - x(i) + x(i + 1) - delta(t, i) * ( K1 + 2 * K2 + 2 * K3 + K4 ) / 6 ;
        end;
    end
    
else
    throw (MException ('ArgumentCheck:IllegalArgument', 'Unsupported method!'));
end

% no inequlity constraints
ci = [];
cij = [];
cih = [];


if ( cejaLength ~= length(ceja) )
    fprintf('Warning: initial size of ceja array: %i differs from end length: %i\n', cejaLength, length(ceja));
end

if ( cehaLength ~= length(ceha) )
    fprintf('Warning: initial size of ceha array: %i differs from end length: %i\n', cehaLength, length(ceha));
end

if ( ~isempty(cej) )
    %  spy(cej);
end

if ( ~isempty(cej) )
    if ( options.checkJacobian && ~isempty(f) )
        optionsCopy = struct(options);
        optionsCopy.checkJacobian = false;
        optionsCopy.checkHessian = false;
        ddeParamEst_constraints_checkCej( cej, f, t, x, n, np, origDelays, xHistory, method, optionsCopy);
    end
end

if ( ~isempty(ceh) )
    if ( options.checkHessian && ~isempty(f) )
        optionsCopy = struct(options);
        optionsCopy.checkJacobian = false;
        optionsCopy.checkHessian = false;
        ddeParamEst_constraints_checkCeh( ceh, f, t, x, n, np, origDelays, xHistory, method, optionsCopy, lm);
    end
end

end