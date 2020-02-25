% dimension
n = [10,5];
fprintf('linear array of %ix%i cells\n',n(1),n(2))
% grid (regular for the moment)
[x,y] = meshgrid(1:n(1),1:n(2));
xx = x'; yy = y';

% # x continuity constraints
nt1 = (n(1)-1)*n(2);
% # y continuity constraints
nt2 = n(1)*(n(2)-1);
% # ground flux conditions
nt3 = n(1);
c_constraints = nt1+nt2+nt3;

% # fluxes per element
factor = 4;
% initialize matrices
B=sparse(factor*prod(n));
D=sparse(prod(n),factor*prod(n));
A=sparse(factor*prod(n));
E=sparse(size(n,2)*prod(n),factor*prod(n));
t1 = sparse(nt1,factor*prod(n));
t2 = sparse(nt2,factor*prod(n));
t3 = sparse(nt3,factor*prod(n));
v0f = zeros(factor*prod(n),1);
% make them full if small
if prod(n)<200, B=full(B);D=full(D);E=full(E);A=full(A);t1=full(t1);t2=full(t2);t3=full(t3); end
% create block matrices
for i=1:prod(n)
    s=(i-1)*factor+1:i*factor; % span of local dofs
    t=(i-1)*size(n,2)+1:i*size(n,2);
    A(s,s)=diag([1,1,2,2]);
    B(s,s)=[-1,0,0,0;0,1,0,0;0,0,-1,0;0,0,0,1];
    D(i,s)=[1,1,1,1];
    E(t,s)=[-.5,.5,0,0;0,0,-.5,.5];
    v0f(s)=[-1,1,-1,1]; % initial wind flux
    % continuity conditions
    [xi,yi]=ind2sub(n,i);    
    if xi > 1
        t1((xi-1)+(yi-1)*(n(1)-1),s(1)) = 1;
    end
    if xi < n(1)
        t1(xi+(yi-1)*(n(1)-1),s(2)) = 1;
    end
    if yi > 1
        t2((yi-1)+(xi-1)*(n(2)-1),s(3)) = 1;
    end
    if yi < n(2)
        t2(yi+(xi-1)*(n(2)-1),s(4)) = 1;
    end
    if yi==1
        t3(xi,s(3)) = 1;
    end
end
C = [t1;t2;t3];
% check that they are the same number computed before
if size(C,1) ~= c_constraints
    error('number of constraints different than indices computed!') 
end
% initial wind flow at the middle of the cells
u0=E*v0f;
% plot initial wind
figure, quiver(xx(:),yy(:),u0(1:2:end),u0(2:2:end)), xlabel('x'), ylabel('y'), title('Initial wind')
% plot initial fluxes
vv0f = B*v0f;
figure, for i=1:prod(n), s=(i-1)*factor+1:i*factor; [xi,yi]=ind2sub(n,i); quiver([xi-.5,xi+.5,xi,xi],[yi,yi,yi-.5,yi+.5],[vv0f(s(1:2))',0,0],[0,0,vv0f(s(3:4))'],.5,'LineWidth',2), hold on, end, xlabel('x'), ylabel('y'), title('Initial wind fluxes')
% need initial Cartesian wind instead of flux
v0=B*v0f; 
% solve mass-consistent using saddle problem
saddle_sparse
% show matrix M if it is small
if prod(n)<10,M,end
% transform fluxes into cartesian winds at the centers
u=E*v;
% plot resulting wind
figure, quiver(xx(:),yy(:),u(1:2:end),u(2:2:end)), xlabel('x'), ylabel('y'), title('Mass-consistent solution')
vv = B*v;
% plot resulting fluxes
figure, for i=1:prod(n), s=(i-1)*factor+1:i*factor; [xi,yi]=ind2sub(n,i); quiver([xi-.5,xi+.5,xi,xi],[yi,yi,yi-.5,yi+.5],[vv(s(1:2))',0,0],[0,0,vv(s(3:4))'],.5,'LineWidth',2), hold on, end, xlabel('x'), ylabel('y'), title('Mass-consistent solution fluxes')
% plot lagrange multiplier p
figure, mesh(xx,yy,reshape(p,n)), xlabel('x'), ylabel('y'), title('Lagrange multiplier p')
% plot lagrange multiplier q
qq = Ct*q;
qqm = (qq(1:4:end-1)+qq(2:4:end)+qq(3:4:end)+qq(4:4:end))/4;
figure, mesh(xx,yy,reshape(qqm,n)), xlabel('x'), ylabel('y'), title('Lagrange multiplier q')
%figure, for i=1:prod(n), s=(i-1)*factor+1:i*factor; [xi,yi]=ind2sub(n,i); mesh([xi-.5,xi+.5;xi,xi],[yi,yi;yi-.5,yi+.5],[qq(s(1:2))';qq(s(3:4))']), hold on, end, xlabel('x'), ylabel('y'), title('Lagrange multiplier q')

