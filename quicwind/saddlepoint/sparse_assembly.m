function [A,D,E,B,C,v0f]=sparse_assembly(X,h)

% sizes and dimensions
d = size(X,2);
n = size(X{1})-1;

% # ground flux conditions
ncg = n(1)*n(2);
% # x continuity constraints
ncx = (n(1)-1)*n(2)*n(3);
% # y continuity constraints
ncy = n(1)*(n(2)-1)*n(3);
% total number of constraints
c_constraints = ncg + ncx + ncy;
if d == 3
    % # z continuity constraints
    ncz = n(1)*n(2)*(n(3)-1);
    c_constraints = c_constraints + ncz;
end

% # fluxes per element
factor = 2*d;

% initialize matrices
A=sparse(factor*prod(n));
D=sparse(prod(n),factor*prod(n));
E=sparse(size(n,2)*prod(n),factor*prod(n));
B=sparse(factor*prod(n));
B2=sparse(factor*prod(n));
cg = sparse(ncg,factor*prod(n));
cx = sparse(ncx,factor*prod(n));
cy = sparse(ncy,factor*prod(n));
if d == 3
    cz = sparse(ncz,factor*prod(n));
end
v0f = zeros(factor*prod(n),1);

% create matrices
for i=1:prod(n)
    [xi,yi,zi]=ind2sub(n,i);
    s=(i-1)*factor+1:i*factor; % span of local dofs
    t=(i-1)*d+1:i*d; % span of dimension coordinates
    % matrix of areas and moduli (u1,u2,v1,v2,w1,w2)
    A(s,s)=diag([1,1,1,1,1e4,1e4]); % v1(1), v2(1), v1(2), v2(2), v1(3), v2(3) 
    % matrix of flux to wind transformation
    B(s,s)=[-1,0,0,0,0,0;   % left face flux to x1 of the wind
            0,1,0,0,0,0;    % right face flux to x1 of the wind
            0,0,-1,0,0,0;   % front face flux to x2 of the wind
            0,0,0,1,0,0;    % back face flux to x2 of the wind
            0,0,0,0,-1,0;   % bottom face flux to x3 of the wind
            0,0,0,0,0,1];   % top face flux to x3 of the wind
    B2(s,s) = cell_B(h(1),h(2),X{3},xi,yi,zi);
    % matrix of divergence free
    D(i,s)=[1,1,1,1,1,1];
    % matrix of resulting winds to winds at the center of the cells
    E(t,s)=[-.5,.5,0,0,0,0;
            0,0,-.5,.5,0,0;
            0,0,0,0,-.5,.5];
    % initial wind flux
    v0f(s)=[-1,1,0,0,1,-1]; 
    % continuity conditions
    if d == 3
        [cx,cy,cz,cg]=c_conditions_3d(n,i,s,cx,cy,cz,cg);
    end
end
C = [cx;cy;cz;cg];
if size(C,1) ~= c_constraints
    error('number of constraints different than indices computed!') 
end
end

function [cx,cy,cz,cg]=c_conditions_3d(n,i,s,cx,cy,cz,cg)
    [xi,yi,zi]=ind2sub(n,i);
    if xi > 1
        cx((xi-1)+(yi-1)*(n(1)-1)+(zi-1)*n(2)*(n(1)-1),s(1)) = 1;
    end
    if xi < n(1)
        cx(xi+(yi-1)*(n(1)-1)+(zi-1)*n(2)*(n(1)-1),s(2)) = 1;
    end
    if yi > 1
        cy((yi-1)+(xi-1)*(n(2)-1)+(zi-1)*n(1)*(n(2)-1),s(3)) = 1;
    end
    if yi < n(2)
        cy(yi+(xi-1)*(n(2)-1)+(zi-1)*n(1)*(n(2)-1),s(4)) = 1;
    end
    if zi > 1
        cz((zi-1)+(xi-1)*(n(3)-1)+(yi-1)*n(1)*(n(3)-1),s(5)) = 1;
    end
    if zi < n(3)
        cz(zi+(xi-1)*(n(3)-1)+(yi-1)*n(1)*(n(3)-1),s(6)) = 1;
    end
    if zi==1
        cg(xi+(yi-1)*n(1),s(5)) = 1;
    end
end