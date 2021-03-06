function [lambda,it,rate,XC]=sparse_solve(K,F,X,params)
n = size(X{1});
nn = size(K,1);
fprintf('sparse_solve: problem size %g mesh %g %g %g\n',nn,n)
rate=0;
XC={};
switch params.solver
    case {'d','direct'}
        disp('sparse direct solver')
        lambda = K\F;
        it=0;
    case {'r','red-black'}
        [lambda,it] = rb_line_gs_solve(K,F,X);
    case {'2','2-level'}
        [lambda,it,rate,XC] = multigrid_solve(K,F,X,params);
    case {'s','schwarz'}
        [lambda,it] = rb_schwarz_solve(K,F,X);
    otherwise
        error('unknown method')
end
end
