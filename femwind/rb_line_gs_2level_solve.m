function [x,it]=rb_line_gs_2level_solve(K,F,X)
% x=rb_line_gs_solve(K,F,X)
coarsening='vertical lines';
coarsening='2';
disp('solver red-black relaxation horizontal, lines vertical')
disp(['coarsening ',coarsening])
n = size(X{1});
nn = size(F,1);
if nn~=prod(n)
    error('rb_line_gs_solve: inconsistent sizes')
end

% settings
tol = 1e-5*big(F)/big(K)

% constant arrays
x = zeros(nn,1);
colx=1:n(3);
onex=ones(1,n(3));

% prolongation
switch coarsening
    case 'vertical lines'
        maxit=500
        nsm = 5
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
    case '2'
        maxit=500
        nsm = 5
        if any(mod(n,2)==0)
            error('number of nodes in each dimension must be odd')
        end
        nc = (n+1)/2-1;
        nnc = prod(nc);  % number of coarse
        P = spalloc(nn,nnc,nnc*27);
        for ic1=1:nc(1)
            for ic2=1:nc(2)
                for ic3=1:nc(3)
                    if1=2*ic1;
                    if2=2*ic2;
                    if3=2*ic3;
                    ixc = sub2ind(nc,ic1,ic2,ic3);
                    for in1=-1:1 
                        for in2=-1:1
                            for in3=-1:1
                                ix=sub2ind(n,if1+in1,if2+in2,if3+in3);
                                P(ix,ixc)=(2-abs(in1))*(2-abs(in2))*(2-abs(in3))/8;
                            end
                        end
                    end
                end
            end
        end
    otherwise
        error(['unknown coarsening ',coarsening])
end
Kc = P'*K*P;
ex = K\F;  % exact solution, for comparison only
for it=1:maxit
    coarse = mod(it,nsm+1)==0;            
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
    semilogy(1:it,res/res(1),'*',1:it,err/err(1),'x',1:it,eer/eer(1),'+'), grid on
    legend('relative 2-norm residual','relative 2-norm error','relative energy norm error')
    title(sprintf('mesh=%g %g %g',n))
    xlabel('iteration')
    drawnow,pause(0.1)
    fprintf('iteration %i residual %g tolerance %g\n',it,res(it),tol)
    if res(it)<tol,
        break
    end
end
end