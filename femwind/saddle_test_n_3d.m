% ordering
% rows: u, v, w.
% columns: low x, high x, low y, high y, low z, high z.

% settings
if ~exist('plot_all','var')
    plot_all = false;
end

% dimension
n = [10,5,5];
h = [1,1,1];
factor = 6;
fprintf('linear array of %ix%ix%i cells\n',n(1),n(2),n(3))

% creating grid
X = regular_mesh(n,h,1);
%thx = h(1)*[0:n(1)]'*ones(1,n(2)+1);
%X = add_terrain_to_mesh(X,thx,'shift');
X = add_terrain_to_mesh(X,'xhill','squash',0.3);
CX = center_mesh(X);
xx = CX{1}; yy = CX{2}; zz = CX{3};

% assembly sparse matrices
[A,D,E,B,C,v0] = sparse_assembly(X,h);

% plot initial wind at the middle of the cells
u0=E*v0;
figure, 
plot_mesh_3d(X), hold on, 
quiver3(xx(:),yy(:),zz(:),u0(1:3:end),u0(2:3:end),u0(3:3:end),'LineWidth',2), xlabel('x'), ylabel('y'), zlabel('z'), title('Initial wind')

if plot_all
    % plot initial wind from each cell
    figure, plot_fluxes_3d(X,B'*v0,h), title('Initial wind fluxes')
end

% solve mass-consistent using saddle problem
saddle_sparse

% transform resulting fluxes into cartesian winds at the centers
u=E*B*v;

% plot resulting wind
figure, 
plot_mesh_3d(X), hold on,  
quiver3(xx(:),yy(:),zz(:),u(1:3:end),u(2:3:end),u(3:3:end),'LineWidth',2), xlabel('x'), ylabel('y'), zlabel('z'), title('Mass-consistent solution')
if plot_all
    % plot resulting fluxes
    figure, plot_fluxes_3d(X,v,h),title('Mass-consistent solution fluxes')
    % plot lagrange multiplier p
    p3 = reshape(p,n);
    figure, for i=1:n(3), mesh(xx(:,:,i),yy(:,:,i),p3(:,:,i)+.05*i); hold on, end, xlabel('x'), ylabel('y'), title('Lagrange multiplier p')
    %{
    figure, for i=1:n(3), qq = Ct*q(:,:,i); qqm = (qq(1:4:end-1)+qq(2:4:end)+qq(3:4:end)+qq(4:4:end))/4; mesh(xx,yy,reshape(qqm,n)), hold on, end, xlabel('x'), ylabel('y'), title('Lagrange multiplier q')
    %}
end