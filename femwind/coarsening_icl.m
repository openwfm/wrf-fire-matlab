function icl=coarsening_icl(X,params)
% [P,X_coarse]=coarsening_2_linear_decide(X,params)
% in:
%   X           grid coordinates
%   params      structure
% out:
%   icl         icl{i} are coarse indices in direction i=1:3

    n = size(X{1});
    nn=prod(n);
    dx=min(X{1}(2:end,1,1)-X{1}(1:end-1,1,1));
    dy=min(X{2}(1,2:end,1)-X{1}(1,1:end-1,1));
    dxy=min(dx,dy);  % horizontal step
    dz = squeeze(X{3}(1,1,2:end)-X{3}(1,1,1:end-1)); % ref z spacing
    % decide on horizontal coarsening factor
    crit=(dz(1)/dxy)/params.a(3);
    if crit > params.minaspect
        hzc=[2,2]; % future proofing if they are different 
    else
        hzc=[1,1];
    end
    hzcavg=sqrt(hzc(1)*hzc(2)); 
    fprintf('horizontal coarsening factor %g %g because weighted height=%g\n',...
        hzc, crit)
    for i=1:2
        icl{i}=1:hzc(i):n(i);
        if icl{i}(end) ~= n(i)
            icl{i}(end+1)=n(i); % last is coarse
        end
        disp(['horizontal ',num2str(i),' coarse layers ',num2str(icl{i})])
    end
    icl3=zeros(1,n(3)); % allocate max
    lcl=1; % last coarse level
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
            icl{3}=icl3(1:i);
            break
        end
    end     
    if nc(3)==0
        error('number of coarse layers is 0')
    end
    disp(['vertical coarse layers ',num2str(icl{3})])
    hg=num2str(X{3}(1,1,icl{3}));
    disp(['heights at corner ',hg])
    hgc=num2str(X{3}(round(n(1)/2),round(n(2)/2),icl{3}));
    disp(['heights at center ',hgc])
    
    nc = [length(icl{1}),length(icl{2}),length(icl{3})];
    disp(['level ',num2str(params.levels),' grid size ',num2str([n,prod(n)]),...
        ' coarse grid size ',num2str([nc,prod(nc)]),' coarsening ratio ',num2str(prod(n)/prod(nc))])
end
