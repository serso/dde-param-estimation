function [ dlmp ] = sqplab_solve_step_ls_iterative(M, me, info, options)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

A = [M, info.ae'; info.ae, zeros(me)];
b = [-info.g; -info.ce];

if ( ~isfield(options, 'dlmp') )
    options.dlmp = [];
end

tol = [];
if ( isfield(options, 'iterativeTol') )
    tol = options.iterativeTol;
end

maxit = [];
if ( isfield(options, 'iterativeMaxit') )
    maxit = options.iterativeMaxit;
end

L = [];
U = [];

%tic;
if ( isfield(options, 'iterativePrecondAlgorithm') )
    iterativePrecondAlgorithmThresh = 0.1;
    if ( isfield(options, 'iterativePrecondAlgorithmThresh') )
        iterativePrecondAlgorithmThresh = options.iterativePrecondAlgorithmThresh;
    end
    
    if ( strcmp(options.iterativePrecondAlgorithm, 'luinc') )
        [L,U] = luinc(A, iterativePrecondAlgorithmThresh);
        if ( strfind(lastwarn,'Incomplete upper triangular factor has') )
            L = [];
            U = [];
        end
    elseif ( strcmp(options.iterativePrecondAlgorithm, 'ilu') )
        setup = [];
        setup.droptol = iterativePrecondAlgorithmThresh;
        [L,U] = ilu(A, setup);
        if ( strfind(lastwarn,'Incomplete upper triangular factor has') )
            L = [];
            U = [];
        end
    elseif ( strcmp(options.iterativePrecondAlgorithm, 'no') )
         % no
    else
        display(options.iterativePrecondAlgorithm);
        throw (MException ('IllegalArgument:NotSupportedMethod', 'Method is not supported'));
    end
end
%toc


%tic;
if ( strcmp(options.stepMethod, 'bicgstab') )
    dlmp = bicgstab(A,b, tol, maxit, L, U, options.dlmp);
elseif ( strcmp(options.stepMethod, 'bicgstabl') )
    dlmp = bicgstabl(A,b, tol, maxit, L, U, options.dlmp);
elseif ( strcmp(options.stepMethod, 'cgs') )
    dlmp = cgs(A,b, tol, maxit, L, U, options.dlmp);
elseif ( strcmp(options.stepMethod, 'gmres') )
    dlmp = gmres(A,b, [], tol, maxit, L, U, options.dlmp);
elseif ( strcmp(options.stepMethod, 'lsqr') )
    dlmp = lsqr(A,b, tol, maxit, L, U, options.dlmp);
elseif ( strcmp(options.stepMethod, 'minres') )
    dlmp = minres(A,b, tol, maxit, L, U, options.dlmp);
elseif ( strcmp(options.stepMethod, 'qmr') )
    dlmp = qmr(A,b, tol, maxit, L, U, options.dlmp);
elseif ( strcmp(options.stepMethod, 'symmlq') )
    dlmp = symmlq(A,b, tol, maxit, L, U, options.dlmp);
elseif ( strcmp(options.stepMethod, 'bicg') )
    dlmp = bicg(A,b, tol, maxit, L, U, options.dlmp);
else 
    display(options.stepMethod);
    throw (MException ('IllegalArgument:NotSupportedMethod', 'Method is not supported'));
end
%toc

if ( isfield(options, 'dlmp') && ~isempty(options.dlmp) )
    % display(norm(abs(dlmp-options.dlmp), inf));
end


end

