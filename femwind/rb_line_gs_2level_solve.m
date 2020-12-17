function [x,it]=rb_line_gs_2level_solve(K,F,X)
% x=rb_line_gs_solve(K,F,X)
disp('solver red-black relaxation horizontal, lines vertical')
n = size(X{1});
nn = size(F,1);
if nn~=prod(n)
    error('rb_line_gs_solve: inconsistent sizes')
end

x = zeros(nn,1);
maxit=250;
colx=1:n(3);
onex=ones(1,n(3));
tol = 1e-5*big(F)/big(K);
ex = K\F;
for it=1:maxit
    for rb1=1:2
        for rb2=1:2
            for i1=rb1:2:n(1)
                for i2=rb2:2:n(2)
                    % solving horizontal location i1 i2 and vertical line
                    ix = sub2ind(n,i1*onex,i2*onex,colx); 
                    x(ix) = x(ix) - K(ix,ix)\(K(:,ix)'*x - F(ix));
                end
            end
        end
    end
    res(it)= norm(K*x-F);
    lambda=zeros(n);
    exact=zeros(n);
    for i=1:nn
        [s1,s2,s3]=ind2sub(n,i);
        lambda(s1,s2,s3)=x(i);
        exact(s1,s2,s3)=ex(i);
    end
    res,err(it)=big(lambda-exact)/big(exact);
    s=7;
    l=squeeze(lambda(:,s,:)-exact(:,s,:));
    xx=squeeze(X{1}(:,s,:));
    yy=squeeze(X{2}(:,s,:));
    zz=squeeze(X{3}(:,s,:));
    figure(3);
    mesh(xx,zz,l)
    t=sprintf('error slice %g y=%g it=%g rel res=%g rel err=%g',...
        s,yy(1),it,res(it)/big(F),err(it)/big(exact));
    title(t)
    figure(4)
    semilogy(1:it,res,'*',1:it,err,'x'), grid on
    legend('residual','error')
    title(sprintf('mesh=%g %g %g',n))
    xlabel('iteration')
    drawnow,pause(0.1)
    fprintf('iteration %g residual %g tolerance %g\n',it,res,tol)
    if res(it)<tol,
        break
    end
end
end