function [P,X_coarse]=coarsening_2_linear(X,params)
% [P,X_coarse]=coarsening_2_linear(X,params)
% in:
%   X           grid coordinates
%   params      structure
% out:
%   P           prolongation matrix, sparse
%   X_coarse    coarse grid coordinates
    n = size(X{1});
    nn=prod(n);
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
    fprintf('horizontal coarsening factor %g %g because weighted height=%g\n',...
        hzc1, hzc2, crit)
    icl1=1:hzc1:n(1);
    if icl1(end) ~= n(1)
        icl1(end+1)=n(1); % last is coarse
    end
    disp(['horizontal 1 coarse layers ',num2str(icl1)])
    icl2=1:hzc2:n(2);
    if icl2(end) ~= n(2)
        icl2(end+1)=n(2); % last is coarse
    end
    disp(['horizontal 2 coarse layers ',num2str(icl2)])
    % decide on vertical coarse levels
    lcl=1; % last coarse level
    icl3=zeros(1,n(3)); % allocate max
    icl3(1)=lcl;
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
            icl3(i+1)=lcl;
        else % at the top already
            nc(3)=i;
            icl3=icl3(1:i);
            break
        end
    end     
    if nc(3)==0
        error('bad number of coarse layers')
    end
    disp(['vertical coarse layers ',num2str(icl3)])
    hg=num2str(X{3}(1,1,icl3));
    disp(['heights at corner ',hg])
    hgc=num2str(X{3}(round(n(1)/2),round(n(2)/2),icl3));
    disp(['heights at center ',hgc])
    
    nc = [length(icl1),length(icl2),length(icl3)];
    disp(['level ',num2str(params.levels),' grid size ',num2str([n,prod(n)]),...
        ' coarse grid size ',num2str([nc,prod(nc)]),' coarsening ratio ',num2str(prod(n)/prod(nc))])
   
    
    disp('building the prolongation')
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
        if3=icl3(ic3);         % the fine number of the coarse layer
        [ifs3,ife3]=ifse(ic3,icl3,nc(3)); % get start and end of support
        fprintf('vertical coarse layer %g at %g contributes to layers %g : %g\n',ic3,if3,ifs3,ife3)
        for ic1=1:nc(1)        % horizontal loops over coarse points
            if1=icl1(ic1);  % fine grid index of the coarse point
            [ifs1,ife1]=ifse(ic1,icl1,nc(1)); % get start and end of support
            % fprintf('coarse x1 %g at %g contributes to %g : %g\n',ic1,if1,ifs1,ife1)
            for ic2=1:nc(2)          
                if2=icl2(ic2);  % fine grid index of the coarse point
                [ifs2,ife2]=ifse(ic2,icl2,nc(2)); % get start and end of support
                % fprintf('coarse x2 %g at %g contributes to %g : %g\n',ic2,if2,ifs2,ife2)
                for l=1:3      % copy coordinates 
                    X_coarse{l}(ic1,ic2,ic3)=X{l}(if1,if2,if3);
                end
                ixc = sub2ind(nc,ic1,ic2,ic3); % index into coarse matrix
                % loop over fine points coarse point ic1 ic2 ic3
                % contributes to
                % coarse point ic1 ic2 ic3 is if1 if2 if3 on the fine grid
                % interpolating from between to (i1 i2 i3) from if1 if2 if3
                % and the next coarse
                for i1=ifs1:ife1
                    for i2=ifs2:ife2
                        for i3=ifs3:ife3
                            ix=sub2ind(n,i1,i2,i3);
                            if params.P_by_x
                                % should adapt from X 
                                if i1>if1
                                    q1=(X{1}(i1,i2,i3)-X{1}(ife1+1,i2,i3))/(X{1}(if1,i2,i3)-X{1}(ife1+1,i2,i3));
                                elseif i1<if1
                                    q1=(X{1}(i1,i2,i3)-X{1}(ifs1-1,i2,i3))/(X{1}(if1,i2,i3)-X{1}(ifs1-1,i2,i3));
                                else
                                    q1=1;
                                end
                                if i2>if2
                                    q2=(X{2}(i1,i2,i3)-X{2}(i1,ife2+1,i3))/(X{2}(i1,if2,i3)-X{2}(i1,ife2+1,i3));
                                elseif i2<if2
                                    q2=(X{2}(i1,i2,i3)-X{2}(i1,ifs2-1,i3))/(X{2}(i1,if2,i3)-X{2}(i1,ifs2-1,i3));
                                else
                                    q2=1;
                                end
                                if i3>if3
                                    q3=(X{3}(i1,i2,i3)-X{3}(i1,i2,ife3+1))/(X{3}(i1,i2,if3)-X{3}(i1,i2,ife3+1));
                                elseif i3<if3
                                    q3=(X{3}(i1,i2,i3)-X{3}(i1,i2,ifs3-1))/(X{3}(i1,i2,if3)-X{3}(i1,i2,ifs3-1));
                                else
                                    q3=1;
                                end
                                val=q1*q2*q3;
                            else
                                val=(1-0.5*abs(i1-if1))...
                                    *(1-0.5*abs(i2-if2))...
                                    *(1-0.5*abs(i3-if3)); % horizontal
                            end
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
    if params.graphics>=2
        check_P(P,X,X_coarse)
    end
end


function [ifs,ife]=ifse(ic,icl,ncn) 
    if ic>1
        ifs=icl(ic-1)+1; % from above previous coarse     
    else
        ifs=icl(ic);     % itself, there is no previous fine layer
    end
    if ic<ncn
        ife=icl(ic+1)-1; % up to under next layer
    else
        ife=icl(ic);     % itself, there is no next layer
    end
end
