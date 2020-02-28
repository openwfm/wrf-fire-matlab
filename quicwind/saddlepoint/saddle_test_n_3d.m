% ordering
% rows: u, v, w.
% columns: low x, high x, low y, high y, low z, high z.

% dimension
n = [11,10,5];
h = [1,1,1];
vstretch = 1.2;
factor = 6;
fprintf('linear array of %ix%ix%i cells\n',n(1),n(2),n(3))

% creating grid
X = regular_mesh(n,h,vstretch);
CX = cell(1,3);
for k = 1:3
    CX{k} = (X{k}(1:end-1,1:end-1,1:end-1)+X{k}(2:end,1:end-1,1:end-1)+X{k}(1:end-1,2:end,1:end-1)+X{k}(1:end-1,1:end-1,2:end)+X{k}(2:end,2:end,1:end-1)+X{k}(2:end,1:end-1,2:end)+X{k}(1:end-1,2:end,2:end)+X{k}(2:end,2:end,2:end))/8;
end
xx = CX{1}; yy = CX{2}; zz = CX{3};

% assembly sparse matrices
[A,D,E,B,C,v0f] = sparse_assembly(X);

% initial Cartesian wind
v0=B*v0f; 

% plot initial wind from each cell
figure, for i=1:prod(n), s=(i-1)*factor+1:i*factor; [xi,yi,zi]=ind2sub(n,i); quiver3([xi-.5,xi+.5,xi,xi,xi,xi],[yi,yi,yi-.5,yi+.5,yi,yi],[zi,zi,zi,zi,zi-.5,zi+.5],[v0(s(1:2))',0,0,0,0],[0,0,v0(s(3:4))',0,0],[0,0,0,0,v0(s(5:6))'],0,'LineWidth',2), hold on, end, xlabel('x'), ylabel('y'), zlabel('z'), title('Initial wind fluxes')
% plot initial wind at the middle of the cells
u0=E*v0f;
figure, quiver3(xx(:),yy(:),zz(:),u0(1:3:end),u0(2:3:end),u0(3:3:end)), xlabel('x'), ylabel('y'), zlabel('z'), title('Initial wind')

% solve mass-consistent using saddle problem
saddle_sparse

% transform resulting fluxes into cartesian winds at the centers
u=E*v;

% plot resulting wind
figure, quiver3(xx(:),yy(:),zz(:),u(1:3:end),u(2:3:end),u(3:3:end)), xlabel('x'), ylabel('y'), zlabel('z'), title('Mass-consistent solution')
vv = B*v;
% plot resulting fluxes
figure, for i=1:prod(n), s=(i-1)*factor+1:i*factor; [xi,yi,zi]=ind2sub(n,i); quiver3([xi-.5,xi+.5,xi,xi,xi,xi],[yi,yi,yi-.5,yi+.5,yi,yi],[zi,zi,zi,zi,zi-.5,zi+.5],[vv(s(1:2))',0,0,0,0],[0,0,vv(s(3:4))',0,0],[0,0,0,0,vv(s(5:6))'],0,'LineWidth',2), hold on, end, xlabel('x'), ylabel('y'), title('Mass-consistent solution fluxes')
% plot lagrange multiplier p
p3 = reshape(p,n);
figure, for i=1:n(3), mesh(xx(:,:,i),yy(:,:,i),p3(:,:,i)+2*i); hold on, end, xlabel('x'), ylabel('y'), title('Lagrange multiplier p')
