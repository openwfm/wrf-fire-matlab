% ordering
% rows: u, v, w.
% columns: low x, high x, low y, high y, low z, high z.

% settings
plot_all = true;

% dimension
n = [2,2,2];
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

% plot initial wind at the middle of the cells
u0=E*v0f;
figure, 
plot_mesh_3d(X), hold on, 
quiver3(xx(:),yy(:),zz(:),u0(1:3:end),u0(2:3:end),u0(3:3:end)), xlabel('x'), ylabel('y'), zlabel('z'), title('Initial wind')

if plot_all
    % plot initial wind from each cell
    figure, plot_fluxes_3d(X,v0), title('Initial wind fluxes')
end

% solve mass-consistent using saddle problem
saddle_sparse

% transform resulting fluxes into cartesian winds at the centers
u=E*v;

% plot resulting wind
figure, 
plot_mesh_3d(X), hold on,  
quiver3(xx(:),yy(:),zz(:),u(1:3:end),u(2:3:end),u(3:3:end)), xlabel('x'), ylabel('y'), zlabel('z'), title('Mass-consistent solution')
if plot_all
    % wind direction fluxes
    vv = B*v;
    % plot resulting fluxes
    figure, plot_fluxes_3d(X,vv),title('Mass-consistent solution fluxes')
    % plot lagrange multiplier p
    p3 = reshape(p,n);
    figure, for i=1:n(3), mesh(xx(:,:,i),yy(:,:,i),p3(:,:,i)+100*i); hold on, end, xlabel('x'), ylabel('y'), title('Lagrange multiplier p')
end