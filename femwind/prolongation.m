function u=prolongation(uc,hzc,icl3,X,params)
% u=prolongation(uc,icl,X,params)
% Multiply by the prolongation matrix
% In:
%   uc      coarse grid vector
%   icl     cell array size 3, indices of coarse grid in the 3 directions
%   X       grid coordinates
%   params
%
% Out:
%   u      fine grid vector 
  
    n = size(X{1});   % grid size
    nn = prod(n);     % number of grid points
    [icl1,icl2]=hzc2icl(hzc,n);    
    nc = [length(icl1),length(icl2),length(icl3)];
    nnc = prod(nc);  % number of coarse grid points
 
    u = zeros(n);  % allocate
    disp('prolongation')
    for ic3=1:nc(3)           % loop over coarse layers    
        if3=icl3(ic3);         % the fine number of the coarse layer
        [ifs3,ife3]=ifse(ic3,icl3,nc(3)); % get start and end of support
        % fprintf('vertical coarse layer %g at %g contributes to layers %g : %g\n',ic3,if3,ifs3,ife3)
        for ic2=1:nc(2)          
            if2=icl2(ic2);  % fine grid index of the coarse point
            [ifs2,ife2]=ifse(ic2,icl2,nc(2)); % get start and end of support
            for i3=ifs3:ife3
                for i2=ifs2:ife2            
                    for ic1=1:nc(1)        % horizontal loops over coarse points
                        if1=icl1(ic1);  % fine grid index of the coarse point
                        [ifs1,ife1]=ifse(ic1,icl1,nc(1)); % get start and end of support
                        % fprintf('coarse x1 %g at %g contributes to %g : %g\n',ic1,if1,ifs1,ife1)
                        % fprintf('coarse x2 %g at %g contributes to %g : %g\n',ic2,if2,ifs2,ife2)
                        % loop over fine points coarse point ic1 ic2 ic3 contributes to
                        % coarse point ic1 ic2 ic3 is if1 if2 if3 on the fine grid
                        % interpolating from between to (i1 i2 i3) from if1 if2 if3
                        % and the next coarse
                        for i1=ifs1:ife1
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
                                u(i1,i2,i3) = u(i1,i2,i3) + q1*q2*q3*uc(ic1,ic2,ic3);
                        end
                    end
                end
            end
        end
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
