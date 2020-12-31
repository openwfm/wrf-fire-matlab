function check_symmetry(A,s,tol)
% check_sym(A,s,tol)
err_A=norm(A-A','fro');
big_A=norm(A,'fro');
tol_A=big_A*length(A)*tol;
if err_A > tol_A,
    fprintf('max nonsymmetry %g > tolerance %g\n',err_A,tol_A);
    error([s,' must be symmetric']);
end
end
