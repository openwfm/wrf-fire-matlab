function [x,it,rate,X_coarse,P]=rb_line_gs_2level_solve(K,F,X,params)
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
        dx=min(X{1}(2:end,1,1)-X{1}(1:end-1,1,1));
        dy=min(X{2}(1,2:end,1)-X{1}(1,1:end-1,1));
        dxy=min(dx,dy);  % horizontal step
        dz = squeeze(X{3}(1,1,2:end)-X{3}(1,1,1:end-1)); % ref z spacing
        % decide on horizontal coarsening factor
        crit=(dz(1)/dxy)/params.a(3);
        if crit > params.minaspect
            hzc1=2;hzc2=2; % future proofing 
        else
            hzc1=1;hzc2=1;
        end
        hzcavg=sqrt(hzc1*hzc2); 
        nc = ceil((n-1)./[hzc1,hzc2,1])+1;
        fprintf('horizontal coarsening factor %g %g because weighted height=%g\n',...
            hzc1, hzc2, crit)
        % decide on vertical coarse levels
        lcl=1; % last coarse level
        icl=zeros(1,nc(3));
        icl(1)=lcl;
        nc(3)=0;
        for i=1:n(3)
            newlcl=lcl+1; % next coarse level by 1
            if lcl+2 <= n(3) 
                crit = ((dz(lcl)+dz(lcl+1))/(2*dxy*hzcavg/2))/params.a(3);
                if crit < params.maxaspect  
                    newlcl=lcl+2; % next coarse level by 2
                end
            end
            lcl = newlcl;
            if lcl <= n(3)
                icl(i+1)=lcl;
            else % at the top already
                nc(3)=i;
                icl=icl(1:i);
                break
            end
        end     
        if nc(3)==0
            error('bad number of coarse layers')
        end
        disp(['vertical coarse layers ',num2str(icl)])
        hg=num2str(X{3}(1,1,icl));
        disp(['heights at corner ',hg])
        hgc=num2str(X{3}(round(n(1)/2),round(n(2)/2),icl));
        disp(['heights at center ',hgc])
        disp(['coarse grid size ',num2str(nc)])
        disp('building the prolongation matrix')
        nnc = prod(nc);  % number of coarse points
        nncz=nnc*27; % estimate number of nonzeros in coarse matrix
        ia=zeros(nncz,1); % preallocate arrays for P matrix entries
        ja=ia;aa=ia;
        % P = spalloc(nn,nnc,nnc*27);
        k=0;
        for l=1:3
            X_coarse{l}=zeros(nc); % preallocate coarse points coordinates
        end
        for ic3=1:nc(3)           % loop over coarse layers    
            if3=icl(ic3);         % the fine level of the coarse layer
            if ic3>1
                ifs3=icl(ic3-1)+1; % from above previous layer     
            else
                ifs3=icl(ic3);     % there is no previous fine layer
            end
            if ic3<nc(3)
                ife3=icl(ic3+1)-1; % up to under next layer
            else
                ife3=icl(ic3);     % there is no next layer
            end
            fprintf('coarse layer %g at %g contributes to layers %g : %g\n',ic3,if3,ifs3,ife3)
            for ic1=1:nc(1)        % horizontal loops over coarse points
                if1=hzc1*ic1-(hzc1-1);  
                if  hzc1 == 1   
                    ifs1=if1;   % no coarsening 
                    ife1=if1;
                elseif if1 > n(1)  % over high boundary   
                    ifs1=n(1);     % fine mesh indices of the coarse point
                    ife1=n(1);
                    if1 = n(1);
                else
                    ifs1=max(1,if1-1);       % start of the support on the fine mesh
                    ife1=min(n(1),if1+1);    % end of the support on the fine mesh
                end
                % fprintf('coarse x1 %g at %g contributes to %g : %g\n',ic1,if1,ifs1,ife1)
                for ic2=1:nc(2)          
                    if2=hzc2*ic2-(hzc2-1); 
                    if hzc2 == 1     
                       ifs2=if2;      % no coarsening
                       ife2=if2;
                    elseif if2 > n(2) % over high boundary
                        ifs2=n(2);   % fine mesh indices of the coarse point
                        if2 = n(2);
                        ife2=n(2);
                    else
                        ifs2=max(1,if2-1);   % start of the support on the fine mesh
                        ife2=min(n(2),if2+1);% end of the support on the fine mesh
                    end
                    % fprintf('coarse x2 %g at %g contributes to %g : %g\n',ic2,if2,ifs2,ife2)
                    for l=1:3      % copy coordinates 
                        X_coarse{l}(ic1,ic2,ic3)=X{l}(if1,if2,if3);
                    end
                    ixc = sub2ind(nc,ic1,ic2,ic3); % index into coarse matrix
                    % loop over fine points coarse point ic1 ic2 ic3 contributes to
                    for i1=ifs1:ife1
                        for i2=ifs2:ife2
                            for i3=ifs3:ife3
                                % should adapt from X 
                                ix=sub2ind(n,i1,i2,i3);
                                val=(1-0.5*abs(i1-if1))...
                                    *(1-0.5*abs(i2-if2))...
                                    *(1-0.5*abs(i3-if3)); % horizontal
                                k=k+1;
                                ia(k)=ix;
                                ja(k)=ixc;
                                aa(k)=val;
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
                for rb1=1:2
                    for rb2=1:2
                        for i1=rb1:2:n(1)
                            for i2=rb2:2:n(2)
                                % solving horizontal location i1 i2 and vertical line
                                for i3=1:n(3)
                                    ix = sub2ind(n,i1,i2,i3); 
                                    x(ix) = x(ix) - K(ix,ix)\(K(:,ix)'*x - F(ix));
                                end
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
