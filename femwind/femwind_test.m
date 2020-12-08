disp('femwind_test')

% dimensions in elements
sc=1; % mesh scale
n = sc*[10,5,5];
h = [1,1,1]/sc;
fprintf('linear array of %ix%ix%i cells\n',n(1),n(2),n(3))
A = diag([1,1,1]);
lambda = zeros(prod(n+1),1); % placeholder solution

% creating the grid
X = regular_mesh(n,h,1.5^(1/sc));
X = add_terrain_to_mesh(X,'xhill','squash',0.3);
CX = center_mesh(X);

% initial wind at the centers of the elements
u0={ones(n),zeros(n),zeros(n)};
figure(1),clf,hold off 
% plot_mesh_3d(X), hold on, 
quiver3(CX{1},CX{2},CX{3},u0{1},u0{2},u0{3},'LineWidth',2), xlabel('x'), ylabel('y'), zlabel('z'), title('Initial wind')
hold off

% assembly sparse matrices
[K,F,~] = sparse_assembly(A,X,u0,lambda);

% dirichlet boundary conditions
[K,F]=apply_boundary_conditions(K,F,X);

% solve the equations
lambda = sparse_solve(K,F,X,'direct');

% gradient of lambda
[~,~,W] = sparse_assembly(A,X,u0,lambda);

% plot resulting wind
figure(2),clf,hold off
% plot_mesh_3d(X), hold on, 
quiver3(CX{1},CX{2},CX{3},u0{1}+W{1},u0{2}+W{2},u0{3}+W{3},...
    'LineWidth',2), xlabel('x'), ylabel('y'), zlabel('z'), title('Final wind')
hold off
