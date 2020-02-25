% dimension
n = [4,3,2];
fprintf('linear array of %ix%ix%i cells\n',n(1),n(2),n(3))
% grid (regular for the moment)
[x,y,z] = meshgrid(1:n(1),1:n(2),1:n(3));
xx = permute(x,[2 1 3]); yy = permute(y,[2 1 3]); zz = permute(z,[2 1 3]);

% # x continuity constraints
nt1 = (n(1)-1)*n(2)*n(3);
% # y continuity constraints
nt2 = n(1)*(n(2)-1)*n(3);
% # z continuity constraints
nt3 = n(1)*n(2)*(n(3)-1);
% # ground flux conditions
nt4 = n(1)*n(2);
c_constraints = nt1+nt2+nt3+nt4;

% # fluxes per element
factor = 6;
% initialize matrices
B=sparse(factor*prod(n));
D=sparse(prod(n),factor*prod(n));
A=sparse(factor*prod(n));
E=sparse(size(n,2)*prod(n),factor*prod(n));
t1 = sparse(nt1,factor*prod(n));
t2 = sparse(nt2,factor*prod(n));
t3 = sparse(nt3,factor*prod(n));
t4 = sparse(nt4,factor*prod(n));
v0f = zeros(factor*prod(n),1);
% make them full if small
if prod(n)<200, B=full(B);D=full(D);E=full(E); end
% create block matrices
for i=1:prod(n)
    s=(i-1)*factor+1:i*factor; % span of local dofs
    t=(i-1)*size(n,2)+1:i*size(n,2);
    A(s,s)=diag([1,1,1,1,1e4,1e4]);
    B(s,s)=[-1,0,0,0,0,0;0,1,0,0,0,0;0,0,-1,0,0,0;0,0,0,1,0,0;0,0,0,0,-1,0;0,0,0,0,0,1];
    D(i,s)=[1,1,1,1,1,1];
    E(t,s)=[-.5,.5,0,0,0,0;0,0,-.5,.5,0,0;0,0,0,0,-.5,.5];
    v0f(s)=[0,0,0,0,1,-1]; % initial wind flux
    % continuity conditions
    [xi,yi,zi]=ind2sub(n,i);
    if xi > 1
        t1((xi-1)+(yi-1)*(n(1)-1)+(zi-1)*n(2)*(n(1)-1),s(1)) = 1;
    end
    if xi < n(1)
        t1(xi+(yi-1)*(n(1)-1)+(zi-1)*n(2)*(n(1)-1),s(2)) = 1;
    end
    if yi > 1
        t2((yi-1)+(xi-1)*(n(2)-1)+(zi-1)*n(1)*(n(2)-1),s(3)) = 1;
    end
    if yi < n(2)
        t2(yi+(xi-1)*(n(2)-1)+(zi-1)*n(1)*(n(2)-1),s(4)) = 1;
    end
    if zi > 1
        t3((zi-1)+(xi-1)*(n(3)-1)+(yi-1)*n(1)*(n(3)-1),s(5)) = 1;
    end
    if zi < n(3)
        t3(zi+(xi-1)*(n(3)-1)+(yi-1)*n(1)*(n(3)-1),s(6)) = 1;
    end
    if zi==1
        t4(xi+(yi-1)*n(1),s(5)) = 1;
    end
end
C = [t1;t2;t3;t4];
if size(C,1) ~= c_constraints
    error('number of constraints different than indices computed!') 
end
% initial wind flow at the middle of the cells
u0=E*v0f;
% plot initial wind
figure, quiver3(xx(:),yy(:),zz(:),u0(1:3:end),u0(2:3:end),u0(3:3:end)), xlabel('x'), ylabel('y'), zlabel('z'), title('Initial wind')
% plot initial fluxes
vv0 = B*v0f;
figure, for i=1:prod(n), s=(i-1)*factor+1:i*factor; [xi,yi,zi]=ind2sub(n,i); quiver3([xi-.5,xi+.5,xi,xi,xi,xi],[yi,yi,yi-.5,yi+.5,yi,yi],[zi,zi,zi,zi,zi-.5,zi+.5],[vv0(s(1:2))',0,0,0,0],[0,0,vv0(s(3:4))',0,0],[0,0,0,0,vv0(s(5:6))'],.5,'LineWidth',2), hold on, end, xlabel('x'), ylabel('y'), zlabel('z'), title('Initial wind fluxes')
% need initial Cartesian wind instead of flux
v0=B*v0f; 
% solve mass-consistent using saddle problem
saddle_sparse
% show matrix M if it is small
if prod(n)<10,M,end
% transform fluxes into cartesian winds at the centers
u=E*v;
% plot resulting wind
figure, quiver3(xx(:),yy(:),zz(:),u(1:3:end),u(2:3:end),u(3:3:end)), xlabel('x'), ylabel('y'), zlabel('z'), title('Mass-consistent solution')
vv = B*v;
% plot resulting fluxes
figure, for i=1:prod(n), s=(i-1)*factor+1:i*factor; [xi,yi,zi]=ind2sub(n,i); quiver3([xi-.5,xi+.5,xi,xi,xi,xi],[yi,yi,yi-.5,yi+.5,yi,yi],[zi,zi,zi,zi,zi-.5,zi+.5],[vv(s(1:2))',0,0,0,0],[0,0,vv(s(3:4))',0,0],[0,0,0,0,vv(s(5:6))'],.5,'LineWidth',2), hold on, end, xlabel('x'), ylabel('y'), title('Mass-consistent solution fluxes')
