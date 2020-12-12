disp('femwind_test')

% dimensions in elements
sc=1; % mesh scale
n = sc*[20,10,5];
sc2=1;
n(1:2)=n(1:2)*sc2
h = [1,1,1]/sc;
fprintf('linear array of %ix%ix%i cells\n',n(1),n(2),n(3))
A = diag([1,1,1]);
lambda = zeros(prod(n+1),1); % placeholder solution

% creating the grid
X = regular_mesh(n,h,1.5^(1/sc));
X = add_terrain_to_mesh(X,'hill','squash',0.3);
CX = center_mesh(X);

% show mesh
hold off
figure
plot_mesh_3d(X), hold on, 
title('The wind mesh, wind vector in centers, lambda in corners')

% initial wind at the centers of the elements
U0={ones(n),zeros(n),zeros(n)};

% show initial wind
figure 
plot_mesh_3d(X,[1,n(1),1,n(2),1,1]), hold on, 
plot_wind_3d(CX,U0)
hold off
title('Initial wind')

% show initial wind
figure 
plot_mesh_3d(X,[1,n(1),1,n(2),1,1]), hold on, 
plot_wind_3d(CX,U0,1)
hold off
title('Initial wind lowest layer')

% assemble sparse system matrix
[K,F,~] = sparse_assembly(A,X,U0,lambda);

% dirichlet boundary conditions
[K,F]=apply_boundary_conditions(K,F,X);

% solve the equations
lambda = sparse_solve(K,F,X,'d');

% assemble final wind
[~,~,W] = sparse_assembly(A,X,U0,lambda);

% plot resulting wind
figure
plot_mesh_3d(X,[1,n(1),1,n(2),1,1]), hold on, 
plot_wind_3d(CX,W)
hold off
title('Final wind')

figure
plot_mesh_3d(X,[1,n(1),1,n(2),1,1]), hold on, 
plot_wind_3d(CX,W,1)
hold off
title('Final wind lowest layer')

condition_number=scond(K)
n
