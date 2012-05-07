function [ dlmp, p, rp, time ] = sqplab_solve_step_ls( M, me, info, options )
%SQPLAB_SOLVE_STEP_LS Summary of this function goes here
%   Detailed explanation goes here
     timerId = tic;
     
     if ( isfield(options, 'stepMethodIterative') && options.stepMethodIterative )
         dlmp = sqplab_solve_step_ls_iterative(M, me, info, options);
         p = [];
         rp = [];
     else
         [dlmp, p, rp] = sqplab_solve_step_ls_direct(M, me, info, options);
     end
     
     time = toc(timerId);

end

