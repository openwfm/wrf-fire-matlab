function [x,it,rate,XC,P]=rb_line_gs_2level_solve(K,F,X,params)
% x=rb_line_gs_solve(K,F,X)
disp(['coarsening method ',params.coarsening])
disp(['smoothing method ',params.smoothing])

n = size(X{1});
nn = size(F,1);
if nn~=prod(n)
    error('rb_line_gs_solve: inconsistent sizes')
end

% settings
tol = params.restol*norm(F);

% constant arrays
x = zeros(nn,1);
colz=1:n(3);
onez=ones(1,n(3));

% prolongation
switch params.coarsening
    case 'vertical lines'
        nnc = n(1)*n(2);  % number of coarse
        P = spalloc(nn,nnc,nn);
        for i1=1:n(1)
            for i2=1:n(2)
                colw=X{3}(i1,i2,end)-X{3}(i1,i2,:);
                ix = sub2ind(n,i1*onez,i2*onez,colz);
                ixc = sub2ind(n(1:2),i1,i2);
                P(ix,ixc)=squeeze(colw);
            end
        end
    case '2 linear'
        if any(mod(n,2)==0)
            error('number of nodes in each dimension must be odd')
        end
        nc = (n+1)/2-1;
        nnc = prod(nc);  % number of coarse
        nncz=nnc*27;
        ia=zeros(nncz,1);
        ja=ia;aa=ia;
        % P = spalloc(nn,nnc,nnc*27);
        k=0;
        for l=1:3,XC{l}=zeros(nc);end
        for ic1=1:nc(1)
            for ic2=1:nc(2)
                for ic3=1:nc(3)
                    if1=2*ic1-1;
                    if2=2*ic2-1;
                    if3=2*ic3-1;
                    ixc = sub2ind(nc,ic1,ic2,ic3);
                    for l=1:3,XC{l}(ic1,ic2,ic3)=X{l}(if1,if2,if3);end
                    for in1=-1:1 
                        for in2=-1:1
                            for in3=-1:1
                                i1=if1+in1;
                                i2=if2+in2;
                                i3=if3+in3;
                                if i1>0 & i2>0 & i3>0,
                                    ix=sub2ind(n,if1+in1,if2+in2,if3+in3);
                                    val=(2-abs(in1))*(2-abs(in2))*(2-abs(in3))/8;
                                    % P(ix,ixc)=val;
                                    k=k+1;
                                    if k>nncz
                                        error('too many nonzeros')
                                    end
                                    ia(k)=ix;
                                    ja(k)=ixc;
                                    aa(k)=val;
                                end
                            end
                        end
                    end
                end
            end
        end
        P = sparse(ia(1:k),ja(1:k),aa(1:k),nn,nnc);
    otherwise
        error(['unknown coarsening ',params.coarsening])
end
disp('coarse matrix')
Kc = P'*K*P;
if params.exact
    disp('exact solution, for comparison only')
    ex = K\F;  
end
cycles=0;
t_cycle='first cycle not complete yet';
for it=1:params.maxit
    coarse = mod(it,params.nsmooth+1)==0;
    if coarse
        fprintf('iteration %g coarse solve\n',it)
        x = x - P*(Kc\(P'*(K*x-F))); 
        it_type='coarse correction';
    else
        fprintf('iteration %g smoothing by %s\n',it,params.smoothing)
        it_type='smoothing';
        switch params.smoothing
            case {'horizontal planes'}
                % disp('red-black vertical, planes horizontal')
                for rb3=1:2
                    for i3=rb3:2:n(3)
                            [planex,planey]=ndgrid(1:n(1),1:n(2));
                            planez = i3*ones(n(1),n(2));
                            % solving horizontal layer
                            ix = sub2ind(n,planex,planey,planez); 
                            x(ix) = x(ix) - K(ix,ix)\(K(:,ix)'*x - F(ix));
                    end
                end
            case {'vertical lines'}
                % disp('red-black relaxation horizontal, lines vertical')
                for rb1=1:2
                    for rb2=1:2
                        for i1=rb1:2:n(1)
                            for i2=rb2:2:n(2)
                                % solving horizontal location i1 i2 and vertical line
                                ix = sub2ind(n,i1*onez,i2*onez,colz); 
                                x(ix) = x(ix) - K(ix,ix)\(K(:,ix)'*x - F(ix));
                            end
                        end
                    end
                end
            case '3D red-black'
                for rb1=1:2
                    for rb2=1:2
                        for rb3=1:2    
                            for i1=rb1:2:n(1)
                                for i2=rb2:2:n(2)
                                    for i3=rb2:2:n(3)
                                        ix = sub2ind(n,i1,i2,i3); 
                                        x(ix) = x(ix) - K(ix,ix)\(K(:,ix)'*x - F(ix));
                                    end
                                end
                            end
                        end
                    end
                end
            otherwise
                error(['smoothing ',params.smoothing,' unknown'])
        end
    end
    r=F-K*x;  % residual
    if params.exact
        e=x-ex;   % error
    else
        e=[];
    end
    res(it)= norm(r);
    if mod(it,params.nsmooth+1)==params.nsmooth
        cycles=cycles+1;
        rate = (res(it)/norm(F))^(1/cycles);
        t_cycle=sprintf('cycle %g avg rate %g',cycles,rate);
        disp(t_cycle)
    end
    tstring=sprintf('it=%g %s %s %s',it,it_type,t_cycle);
    plot_error_slice(e,r,X,tstring,params)
    figure(14)
    semilogy(1:it,res/norm(F),'*')
    legend('relative 2-norm residual')
    if params.exact
        err(it)=norm(e); % l2 error
        eer(it)=sqrt(e'*K*e); % energy norm error
        hold on, semilogy(1:it,err/err(1),'x',1:it,eer/eer(1),'+'), hold off
        legend('relative 2-norm residual','relative 2-norm error','relative energy norm error')
    end
    grid on
    title([sprintf('mesh=%g %g %g ',n),t_cycle])
    xlabel('iteration')
    drawnow,pause(0.1)
    fprintf('iteration %i residual %g tolerance %g\n',it,res(it),tol)

    if res(it)<tol,
        break
    end
end
end