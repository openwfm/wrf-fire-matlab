function [W,rate]=femwind_solve(A,X,U0,params)
    % assemble sparse system matrix
    nel=size(X{1})-1;
    lambda = zeros(prod(nel+1),1); % placeholder solution
    % F = f_assembly_fortran(A,X,U0,lambda,params);
    [K,F,~] = sparse_assembly(A,X,U0,lambda,params);
    if isfield(params,'femwind_fortran_test') && params.femwind_fortran_test
        F_f=read_array('F_f');
        err_f = big(F(:)-F_f(:))
        if err_f > 1e-4
             error('err_f too large')
        end
    end
    

    % dirichlet boundary conditions
    [K,~]=apply_boundary_conditions(K,[],X);   % to K - once
    [~,F]=apply_boundary_conditions([],F,X);   % to F - every time

    % solve the equations
    % [lambda,it] = sparse_solve(K,F,X,'s');
    [lambda,it,rate,XC] = sparse_solve(K,F,X,params);
    format long
    rate

    % assemble final wind
    [~,~,W] = sparse_assembly(A,X,U0,lambda,params);
end
