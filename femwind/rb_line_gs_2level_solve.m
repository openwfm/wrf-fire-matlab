function [x,it,rate,X_coarse,P]=rb_line_gs_2level_solve(K,F,X,params)
% x=rb_line_gs_solve(K,F,X)
disp(['coarsening method ',params.coarsening])
disp(['smoothing method ',params.smoothing])

n = size(X{1});
nn = size(F,1);
if nn~=prod(n)
    error('rb_line_gs_solve: inconsistent sizes')
end
rate = 0;

% settings
tol = params.restol*norm(F);

% constant arrays
x = zeros(nn,1);
colz=1:n(3);
onez=ones(1,n(3));

disp('building prolongation')
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
        [P,X_coarse]=coarsening_2_linear(X,params);
    otherwise
        error(['unknown coarsening ',params.coarsening])
end
disp('computing coarse matrix')
diary;diary
K_coarse = P'*K*P;
check_nonzeros(params.levels-1,K_coarse,X_coarse,P,K,X);
if params.exact
    disp('exact solution, for comparison only')
    ex = K\F;  
end
cycles=0;
t_cycle='first cycle not complete yet';
for it=1:params.maxit
    diary;diary
    params.it(params.levels)=it;  % where we are, for prints and filenames
    coarse = mod(it,params.nsmooth+1)==0;
    if coarse
        fprintf('iteration %g level %g coarse correction\n',it,params.levels)
        F_coarse = P'*(K*x-F);
        if params.apply_coarse_boundary_conditions
            [K_coarse,F_coarse]=apply_boundary_conditions(K_coarse,F_coarse,X_coarse);
        end
        if params.levels<=2 % next is 1, the coarsest
            x_coarse = K_coarse\F_coarse;
        else  % solve coarse problem recursively
            params_coarse=params;  % copy all params 
            params_coarse.levels=params.levels-1;
            params_coarse.nsmooth=params.nsmooth_coarse;
            params_coarse.maxit=params.maxit_coarse;
            params_coarse.iterations_fig=params.iterations_fig+10;
            params_coarse.res_slice_fig=params.res_slice_fig+10;
            params_coarse.err_slice_fig=params.err_slice_fig+10;
            [x_coarse,~,~,~,~]=rb_line_gs_2level_solve(K_coarse,F_coarse,X_coarse,params_coarse);
        end
        fprintf('coarse solve done, level %g continuting\n',params.levels)
        x = x - P*x_coarse; 
        it_type=sprintf('level %g coarse correction',params.levels);
    else
        fprintf('iteration %g level %g smoothing by %s\n',it,params.levels,params.smoothing)
        it_type=sprintf('level %g smoothing',params.levels);
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
            case {'vertical sweeps'}
                % disp('red-black relaxation horizontal, down to up sweep vertical')
                x = vertical_sweeps(K,F,X,x);
            case '3D red-black'
                for rb1=1:2
                    for rb2=1:2
                        for rb3=1:2    
                            for i1=rb1:2:n(1)
                                for i2=rb2:2:n(2)
                                    for i3=rb3:2:n(3)
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
    r=F-K*x;  % residualm diagnostics only
    if params.exact
        e=x-ex;   % error
    else
        e=[];
    end
    res(it)= norm(r);
    if mod(it,params.nsmooth+1)==params.nsmooth
        cycles=cycles+1;
        rate = (res(it)/norm(F))^(1/cycles);
        t_cycle=sprintf('cycle %g level %g avg rate %g',cycles,params.levels,rate);
        disp(t_cycle)
    end
    tstring=sprintf('it=%g %s %s %s',it,it_type,t_cycle);
    plot_error_slice(e,r,X,tstring,params)
    figure(params.iterations_fig)
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
if params.save_files > 2
    sfile=sprintf('%s_%s_%i.mat',params.save_file_prefix,params.id,params.levels);
    save(sfile,'-v7.3')
end
end
