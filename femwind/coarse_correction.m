function x=coarse_correction(x,F,K,K_coarse,X_coarse,hzc,icl3,X,params)
    n = size(X{1});
    F_coarse3d = restriction(reshape(F-K*x,n),hzc,icl3,X,params);
    F_coarse = F_coarse3d(:);
    if params.apply_coarse_boundary_conditions
        [K_coarse,F_coarse]=apply_boundary_conditions(K_coarse,F_coarse,X_coarse);
    end
    if params.levels<=2 % next is 1, the coarsest
        if params.coarsest_iter==0
            x_coarse = K_coarse\F_coarse;
        else
            x_coarse =zeros(size(F_coarse));
            for it=1:params.coarsest_iter
                x_coarse=smoothing(K_coarse,F_coarse,X_coarse,x_coarse,params);
            end
        end
    else  % solve coarse problem recursively
        params_coarse=params;  % copy all params 
        params_coarse.levels=params.levels-1;
        params_coarse.nsmooth=params.nsmooth_coarse;
        params_coarse.maxit=params.maxit_coarse;
        params_coarse.iterations_fig=params.iterations_fig+10;
        params_coarse.res_slice_fig=params.res_slice_fig+10;
        params_coarse.err_slice_fig=params.err_slice_fig+10;
        [x_coarse,~,~,~]=multigrid_solve(K_coarse,F_coarse,X_coarse,params_coarse);
    end
    fprintf('coarse solve done, level %g continuting\n',params.levels)
    % x = x + P*x_coarse
    x_increment = prolongation(reshape(x_coarse,size(X_coarse{1})),hzc,icl3,X,params);
    x = x + x_increment(:);
end