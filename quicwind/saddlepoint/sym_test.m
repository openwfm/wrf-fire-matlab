%% symbolic settings
% dimensions and symbolic step sizes
n = [5,4,3];
h = [sym('dx'),sym('dy'),sym('dz')];

%% regular_mesh function
val = 1;
zz = sym(zeros(1,n(3)));
zz(1) = 0;
for i=2:n(3)+1
    zz(i) = zz(i-1) + h(3) * val^(i-2);
end
[x,y,z] = ndgrid(h(1)*[0:n(1)],h(2)*[0:n(2)],zz);
X = {x,y,z};

%% add_terrain_to_mesh function
thx = h(1)*[0:n(1)]'*ones(1,n(2)+1);
kmax=size(X{1},3);
XX=X;
for k=1:kmax
    XX{3}(:,:,k)=X{3}(:,:,k)+thx;
end
X=XX;

%% sparse_assembly
% moduli for penalization
moduli = [sym('alpha_1'),sym('alpha_1'),sym('alpha_2')];
% sizes and dimensions
d = size(X,2);
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
A = sym(sparse(factor*prod(n)));
B = sym(sparse(factor*prod(n)));
D = sparse(prod(n),factor*prod(n));
E = sparse(size(n,2)*prod(n),factor*prod(n));
cg = sparse(ncg,factor*prod(n));
cx = sparse(ncx,factor*prod(n));
cy = sparse(ncy,factor*prod(n));
cz = sparse(ncz,factor*prod(n));
v0 = sym(zeros(factor*prod(n),1));
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
    v0(s)=[sym('vx_1'),sym('vx_2'),sym('vy_1'),sym('vy_2'),sym('vz_1'),sym('vz_2')]';
    % continuity conditions
    [cx,cy,cz,cg]=c_conditions_3d(n,xi,yi,zi,s,cx,cy,cz,cg);
end
% continuity operator
C = [cx;cy;cz;cg];
% check number of continuity constraints
if size(C,1) ~= c_constraints
    error('number of constraints different than indices computed!') 
end