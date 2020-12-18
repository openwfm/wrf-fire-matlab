function [x,it]=rb_line_gs_2level_solve(K,F,X)
% x=rb_line_gs_solve(K,F,X)
disp('solver red-black relaxation horizontal, lines vertical')
n = size(X{1});
nn = size(F,1);
if nn~=prod(n)
    error('rb_line_gs_solve: inconsistent sizes')
end

% settings
maxit=100
nsm = 10
tol = 1e-5*big(F)/big(K)

% constant arrays
x = zeros(nn,1);
colx=1:n(3);
onex=ones(1,n(3));

% prolongation
nnc = n(1)*n(2);  % number of coarse
P = spalloc(nn,nnc,nn);
for i1=1:n(1)
    for i2=1:n(2)
        colw=X{3}(i1,i2,end)-X{3}(i1,i2,:);
        ix = sub2ind(n,i1*onex,i2*onex,colx);
        ixc = sub2ind(n(1:2),i1,i2);
        P(ix,ixc)=squeeze(colw);
    end
end
Kc = P'*K*P;
ex = K\F;  % exact solution, for comparison only
for it=1:maxit
    coarse = mod(it,nsm+1),            
    if coarse
        
        x = x - P*(Kc\(P'*(K*x-F))); % coarse solve
    else  % smoothing
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
    end
    res(it)= norm(K*x-F);
    lambda=zeros(n);
    exact=zeros(n);
    for i=1:nn
        [s1,s2,s3]=ind2sub(n,i);
        lambda(s1,s2,s3)=x(i);
        exact(s1,s2,s3)=ex(i);
    end
    err(it)=norm(x-ex); % l2 error
    eer(it)=sqrt((x-ex)'*K*(x-ex)); % energy norm error
    s=7;
    l=squeeze(lambda(:,s,:)-exact(:,s,:));
    xx=squeeze(X{1}(:,s,:));
    yy=squeeze(X{2}(:,s,:));
    zz=squeeze(X{3}(:,s,:));
    figure(13);
    mesh(xx,zz,l)
    t=sprintf('error slice %g y=%g it=%g rel res=%g rel err=%g eer=%g',...
        s,yy(1),it,res(it)/res(1),err(it)/err(1),eer(it)/eer(1));
    xlabel('horizontal')
    ylabel('vertical')
    title(t)
    figure(14)
    semilogy(1:it,res,'*',1:it,err,'x',1:it,eer,'+'), grid on
    legend('relative 2-norm residual','relative 2-norm error','relative energy norm error')
    title(sprintf('mesh=%g %g %g',n))
    xlabel('iteration')
    drawnow,pause(0.1)
    fprintf('iteration %g residual %g tolerance %g\n',it,res,tol)
    if res(it)<tol,
        break
    end
end
end