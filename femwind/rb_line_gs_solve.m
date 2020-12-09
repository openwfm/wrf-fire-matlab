function x=rb_line_gs_solve(K,F,X)
% x=rb_line_gs_solve(K,F,X)
disp('solver red-black relaxation horizontal, lines vertical')
n = size(X{1});
nn = size(F,1);
if nn~=prod(n)
    error('rb_line_gs_solve: inconsistent sizes')
end

x = zeros(nn,1);
maxit=1000;
colx=1:n(3);
onex=ones(1,n(3));
tol = 1e-5*big(F)/big(K);
for it=1:maxit
    for rb1=1:2
        for rb2=1:2
            for i1=rb1:2:n(1)
                for i2=rb2:2:n(2)
                    % solving horizontal location i1 i2 and vertical line
                    ix = sub2ind(n,i1*onex,i2*onex,colx); 
                    x(ix) = x(ix) - K(ix,ix)\(K(ix,:)*x - F(ix));
                end
            end
        end
    end
    res= norm(K*x-F);
    fprintf('iteration %g residual %g tolerance %g\n',it,res,tol)
    if res<tol,
        break
    end
end