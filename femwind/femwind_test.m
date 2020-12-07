disp('femwind_test')

% dimensions in elements
n = [10,5,5];
h = [1,1,1];
fprintf('linear array of %ix%ix%i cells\n',n(1),n(2),n(3))
A = diag([1,1,2]);

% creating grid
X = regular_mesh(n,h,1.5);
X = add_terrain_to_mesh(X,'xhill','squash',0.3);
CX = center_mesh(X);

% initial wind at the centers of the elements
u0={ones(n),zeros(n),zeros(n)};
figure, 
plot_mesh_3d(X), hold on, 
quiver3(CX{1},CX{2},CX{3},u0{1},u0{2},u0{3},'LineWidth',2), xlabel('x'), ylabel('y'), zlabel('z'), title('Initial wind')

% assembly sparse matrices
[K,F] = sparse_assembly(A,X,u0);

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