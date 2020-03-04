function [A,D,E,B,C,v0]=sparse_assembly(X,h)

% sizes and dimensions
d = size(X,2);
n = size(X{1})-1;

% # ground flux conditions
ncg = n(1)*n(2);
% # x continuity constraints
ncx = (n(1)-1)*n(2)*n(3);
% # y continuity constraints
ncy = n(1)*(n(2)-1)*n(3);
% # z continuity constraints
ncz = n(1)*n(2)*(n(3)-1);
% total number of constraints
c_constraints = ncg + ncx + ncy + ncz;

% # fluxes per element
factor = 2*d;

% initialize matrices
A=sparse(factor*prod(n));
D=sparse(prod(n),factor*prod(n));
E=sparse(size(n,2)*prod(n),factor*prod(n));
B=sparse(factor*prod(n));
cg = sparse(ncg,factor*prod(n));
cx = sparse(ncx,factor*prod(n));
cy = sparse(ncy,factor*prod(n));
cz = sparse(ncz,factor*prod(n));
v0 = zeros(factor*prod(n),1);

% moduli for penalization
moduli = [1e3,1e3,1];

% create matrices
for i=1:prod(n)
    [xi,yi,zi]=ind2sub(n,i);
    s=(i-1)*factor+1:i*factor; % span of local dofs
    t=(i-1)*d+1:i*d; % span of dimension coordinates
    % matrix A of areas and moduli (u1,u2,v1,v2,w1,w2), and
    % matrix B of flux to wind transformation
    [A(s,s),B(s,s)] = cell_A_B(h(1),h(2),X{3},xi,yi,zi,moduli);
    % matrix of divergence operator
    D(i,s)=[1,1,1,1,1,1];
    % matrix of resulting winds to winds at the center of the cells
    E(t,s)=[.5,.5,0,0,0,0;
            0,0,.5,.5,0,0;
            0,0,0,0,.5,.5];
    % initial wind
    v0(s)=[1,1,0,0,0,0]';
    % continuity conditions
    [cx,cy,cz,cg]=c_conditions_3d(n,i,s,cx,cy,cz,cg);
end
% continuity operator
C = [cx;cy;cz;cg];
% check number of continuity constraints
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