function C=midpoints2corners(M)
% extrapolation formula
extra = @(a,b) 2*a-b;

n=size(M);
d=length(n);
C=zeros(n+1);
if d==1
    % inner-corners
    C(2:end-1)=avg1d(M);
    % corners
    C(1)=extra(C(2),C(3));
    C(end)=extra(C(end-1),C(end-2));
elseif d==2
    % inner-corners
    C(2:end-1,2:end-1)=avg2d(M);
    % edge-corners
    C(1,2:end-1)=extra(C(2,2:end-1),C(3,2:end-1));
    C(end,2:end-1)=extra(C(end-1,2:end-1),C(end-2,2:end-1));
    C(2:end-1,1)=extra(C(2:end-1,2),C(2:end-1,3));
    C(2:end-1,end)=extra(C(2:end-1,end-1),C(2:end-1,end-2));
    % corners
    C(1,1)=extra(C(2,2),C(3,3)); 
    C(1,end)=extra(C(2,end-1),C(3,end-2));
    C(end,1)=extra(C(end-1,2),C(end-2,3));
    C(end,end)=extra(C(end-1,end-1),C(end-2,end-2));
elseif d==3
    % inner-corners
    C(2:end-1,2:end-1,2:end-1)=avg3d(M);
    % face-corners
    C(1,2:end-1,2:end-1)=extra(C(2,2:end-1,2:end-1),C(3,2:end-1,2:end-1));
    C(end,2:end-1,2:end-1)=extra(C(end-1,2:end-1,2:end-1),C(end-2,2:end-1,2:end-1));
    C(2:end-1,1,2:end-1)=extra(C(2:end-1,2,2:end-1),C(2:end-1,3,2:end-1));
    C(2:end-1,end,2:end-1)=extra(C(2:end-1,end-1,2:end-1),C(2:end-1,end-2,2:end-1));
    C(2:end-1,2:end-1,1)=extra(C(2:end-1,2:end-1,2),C(2:end-1,2:end-1,3));
    C(2:end-1,2:end-1,end)=extra(C(2:end-1,2:end-1,end-1),C(2:end-1,2:end-1,end-2));
    % edge-corners
    C(1,1,2:end-1)=extra(C(2,2,2:end-1),C(3,3,2:end-1));
    C(end,1,2:end-1)=extra(C(end-1,2,2:end-1),C(end-2,3,2:end-1));
    C(1,end,2:end-1)=extra(C(2,end-1,2:end-1),C(3,end-2,2:end-1));
    C(end,end,2:end-1)=extra(C(end-1,end-1,2:end-1),C(end-2,end-2,2:end-1));
    C(1,2:end-1,1)=extra(C(2,2:end-1,2),C(3,2:end-1,3));
    C(end,2:end-1,1)=extra(C(end-1,2:end-1,2),C(end-2,2:end-1,3));
    C(1,2:end-1,end)=extra(C(2,2:end-1,end-1),C(3,2:end-1,end-2));
    C(end,2:end-1,end)=extra(C(end-1,2:end-1,end-1),C(end-2,2:end-1,end-2));
    C(2:end-1,1,1)=extra(C(2:end-1,2,2),C(2:end-1,3,3));
    C(2:end-1,end,1)=extra(C(2:end-1,end-1,2),C(2:end-1,end-2,3));
    C(2:end-1,1,end)=extra(C(2:end-1,2,end-1),C(2:end-1,3,end-2));
    C(2:end-1,end,end)=extra(C(2:end-1,end-1,end-1),C(2:end-1,end-2,end-2));
    % corners
    C(1,1,1)=extra(C(2,2,2),C(3,3,3));
    C(end,1,1)=extra(C(end-1,2,2),C(end-2,3,3));
    C(1,end,1)=extra(C(2,end-1,2),C(3,end-2,3));
    C(1,1,end)=extra(C(2,2,end-1),C(3,3,end-2));
    C(1,end,end)=extra(C(2,end-1,end-1),C(3,end-2,end-2));
    C(end,1,end)=extra(C(end-1,2,end-1),C(end-2,3,end-2));
    C(end,end,1)=extra(C(end-1,end-1,2),C(end-2,end-2,3));
    C(end,end,end)=extra(C(end-1,end-1,end-1),C(end-2,end-2,end-2));
else
    error(['Not implemented number of dimensions ',d]);
end


function A=avg1d(D)
    A=( D(1:end-1) + D(2:end) ) / 2;
end

function A=avg2d(D)
    A=( D(1:end-1,1:end-1) + D(1:end-1,2:end) + ...
        D(2:end,1:end-1) + D(2:end,2:end) ) / 4;
end

function A=avg3d(D)
    A=( D(1:end-1,1:end-1,1:end-1) + D(1:end-1,1:end-1,2:end) + ...
        D(1:end-1,2:end,1:end-1) + D(1:end-1,2:end,2:end) + ...
        D(2:end,1:end-1,1:end-1) + D(2:end,1:end-1,1:end-1) + ...
        D(1:end-1,1:end-1,1:end-1) + D(1:end-1,1:end-1,1:end-1) ) / 8;
end
end
