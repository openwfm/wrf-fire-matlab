function [hzc,icl3]=coarsening_icl(dx,dy,dz,params)
% [[hzc,icl3]=coarsening_icldx,dy,dz,params)
% in:
%   dx,dy       mesh spacings, scalars
%   dz          vertical element sizes, vector
%   params      structure
% out:
%   hzc         horizontal coarsening factors in directions 1 and 2
%   icl3        coarse indices in direction 3
% 
    if ~isvector(dz),
          error('dz must be a vector')
    end
    dz = dz(:)';  % make sure dz is a row
    disp(['coarsening_icl: input dx,dy,dz=',num2str([dx,dy,dz])])
    dxy=min(dx,dy);  % horizontal step
    n3 = length(dz)+1;
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
    icl3=zeros(1,n3); % allocate max
    lcl=1; % last coarse level
    icl3(1)=lcl;
    nc3=0;
    for i=1:n3
        newlcl=lcl+1; % next coarse level by 1
        if lcl+2 <= n3 
            crit = ((dz(lcl)+dz(lcl+1))/(2*dxy*hzcavg/2))/params.a(3);
            if crit < params.maxaspect  
                newlcl=lcl+2; % next coarse level by 2
            end
        end
        lcl = newlcl;
        if lcl <= n3
            icl3(i+1)=lcl;
        else % at the top already
            nc3=i;
            icl3 = icl3(1:i);
            break
        end
    end     
    if nc3==0
        error('number of coarse layers is 0')
    end
    disp(['vertical coarse layers ',num2str(icl3)])
    hg=[0,cumsum(dz)];
    hgc=hg(icl3);
    disp(['heights above terrain ',num2str(hg)])
    disp(['coarse heights above terrain ',num2str(hgc)])
end
