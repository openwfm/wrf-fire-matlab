function [x,it,rate,X_coarse]=multigrid_solve(K,F,X,params)
% x=multigrid_solve(K,F,X)

n = size(X{1});
fprintf('multigrid level %i grid size %i %i %i rhs %g\n',params.levels,n,big(F))

disp(['coarsening method ',params.coarsening])
disp(['smoothing method ',params.smoothing])

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
        % [P,X_coarse]=coarsening_2_linear(X,params);
        dx=(X{1}(2,1,1)-X{1}(1,1,1));
        dy=(X{2}(1,2,1)-X{1}(1,1,1));
        dz = squeeze(X{3}(1,1,2:end)-X{3}(1,1,1:end-1)); % ref z spacing
        fprintf('layers: '),disp(X{3}(1,1,:)),fprintf('dz='),disp(dz)
        [hzc,icl3]=coarsening_icl_fortran(dx,dy,dz,params);
        X_coarse=coarsening_X(hzc,icl3,X,params);
        nnc=prod(size(X_coarse{1}));
        % P=coarsening_P(icl,X,params);
        prol_rest_err(hzc,icl3,X,params);
    otherwise
        error(['unknown coarsening ',params.coarsening])
end
fprintf('coarsening level variables %i coarse %i ratio %g\n',params.levels,nn,nnc,nn/nnc);  
disp('computing coarse matrix')
diary;diary
switch params.coarse_K
    case {1,'variational'}
        disp('coarse K is variational')
        P=coarsening_P(icl,X,params);
        K_coarse = P'*K*P;
        check_nonzeros(params.levels-1,K_coarse,X_coarse,P,K,X);
    case {2,'assembly'}
        disp('coarse K by assembly')
        % assemble sparse system matrix
        [K_coarse,~,~] = sparse_assembly(diag(params.a),X_coarse,[],[],params);
        % dirichlet boundary conditions
        [K_coarse,~]=apply_boundary_conditions(K_coarse,[],X_coarse);
    otherwise
        disp(params.coarse_P)
        error('unknown coarse_P method')
end
if params.exact
    disp('exact solution, for comparison only')
    ex = K\F;  
end
cycles=0;
t_cycle='first cycle not complete yet';
if params.levels<=1 || nn == nnc % coarsest level
    if params.coarsest_iter==0   % direct 
        fprintf('multigrid coarsest level %i solving directly\n',params.levels)
        x = K\F;
    else
        x =zeros(size(F));       % iterative
        fprintf('multigrid coarsest level %i solving by %i iterations\n',params.levels,params.coarsest_iter)
        for it=1:params.coarsest_iter
            x=smoothing(K,F,X,x,params);
        end
        res=big(K*x-F);relres=res/big(F);rate=relres^(1/params.coarsest_iter);
        fprintf('multigrid coarsest residual %g relative %g rate %g\n',res,relres,rate);
    end
    return
end
for it=1:params.maxit
    diary;diary
    params.it(params.levels)=it;  % where we are, for prints and filenames
    coarse = mod(it,params.nsmooth+1)==0;
    if coarse
        fprintf('iteration %g level %g coarse correction\n',it,params.levels)
        x=coarse_correction(x,F,K,K_coarse,X_coarse,hzc,icl3,X,params);
        it_type=sprintf('level %g coarse correction',params.levels);
    else
        it_type='smoothing';
        fprintf('iteration %g level %g smoothing by %s\n',it,params.levels,params.smoothing)
        x=smoothing(K,F,X,x,params);
    end
    r=F-K*x;  % residual for diagnostics only
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
    if params.graphics > -1
        tstring=sprintf('it=%g %s %s %s',it,it_type,t_cycle);
        plot_error_slice(e,r,X,tstring,params)
    end
    if params.graphics > -2
        figure(params.iterations_fig)
        semilogy(1:it,res/norm(F),'*')
        legend('relative 2-norm residual')
        grid on
        title([sprintf('mesh=%g %g %g ',n),t_cycle])
        xlabel('iteration')
        drawnow,pause(0.1)
    end
    if params.exact
        err(it)=norm(e); % l2 error
        eer(it)=sqrt(e'*K*e); % energy norm error
        hold on, semilogy(1:it,err/err(1),'x',1:it,eer/eer(1),'+'), hold off
        legend('relative 2-norm residual','relative 2-norm error','relative energy norm error')
    end
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


