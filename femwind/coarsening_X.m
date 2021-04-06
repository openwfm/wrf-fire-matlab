function X_coarse=coarsening_X(icl,X,params)
% [P,X_coarse]=coarsening_2_linear_decide(X,params)
% in:
%   icl         icl{i} are coarse indices in direction i=1:3
%   X           grid coordinates
%   params      structure
% out:
%   X_coarse    coarse grid coordinates

    n = size(X{1});
    nn=prod(n);
    icl1=icl{1};icl2=icl{2};icl3=icl{3};  % unwrap cell array
    nc = [length(icl1),length(icl2),length(icl3)];
    for l=1:3
        X_coarse{l}=zeros(nc); % preallocate coarse points coordinates
    end
    for ic3=1:nc(3)           % loop over coarse layers    
        if3=icl3(ic3);         % the fine number of the coarse layer
        for ic1=1:nc(1)        % horizontal loops over coarse points
            if1=icl1(ic1);  % fine grid index of the coarse point
            for ic2=1:nc(2)          
                if2=icl2(ic2);  % fine grid index of the coarse point
                for l=1:3      % copy coordinates 
                    X_coarse{l}(ic1,ic2,ic3)=X{l}(if1,if2,if3);
                end
            end
        end
    end
end
