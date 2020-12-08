function lambda=sparse_solve(K,F,X,method)
n = size(X{1});
nn = size(K,1);
fprintf('sparse_solve: problem size %g mesh %g %g %g\n',nn,n)
switch method
    case {'d','direct'}
        disp('sparse direct solver')
        lambda = K\F;
    case {'r','red-black'}
        lambda=rb_line_gs_solve(K,F,X);
    otherwise
        error('unknown method')
end
end
