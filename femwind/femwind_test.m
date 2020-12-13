disp('femwind_test')

do_plot=1

% dimensions in elements
sc=1; % mesh scale
n = sc*[20,10,5];
sc2=1;
n(1:2)=n(1:2)*sc2
h = [1,1,1]/sc;
fprintf('linear array of %ix%ix%i cells\n',n(1),n(2),n(3))
da=[1 1 1]
string_diag_A=sprintf('%g %g %g',da);
A = diag(da);
lambda = zeros(prod(n+1),1); % placeholder solution

% creating the grid
X = regular_mesh(n,h,1.5^(1/sc));
X = add_terrain_to_mesh(X,'hill','squash',0.3);
CX = center_mesh(X);

% initial wind at the centers of the elements
U0={ones(n),zeros(n),zeros(n)};

if do_plot

    % show mesh
    hold off
    figure
    plot_mesh_3d(X), hold on, 
    title('The wind mesh, wind vector in centers, lambda in corners')

    % show initial wind
    figure 
    plot_mesh_3d(X,[1,n(1),1,n(2),1,1]), hold on, 
    plot_wind_3d(CX,U0)
    hold off
    axis equal
    title('Initial wind')

    % show initial wind
    figure 
    plot_mesh_3d(X,[1,n(1),1,n(2),1,1]), hold on, 
    plot_wind_3d(CX,U0,1)
    hold off
    axis equal
    title('Initial wind lowest layer')
end

% assemble sparse system matrix
[K,F,~] = sparse_assembly(A,X,U0,lambda);

% dirichlet boundary conditions
[K,F]=apply_boundary_conditions(K,F,X);

% solve the equations
lambda = sparse_solve(K,F,X,'s');

% assemble final wind
[~,~,W] = sparse_assembly(A,X,U0,lambda);

if do_plot
    % plot resulting wind
    figure
    plot_mesh_3d(X,[1,n(1),1,n(2),1,1]), hold on, 
    plot_wind_3d(CX,W)
    hold off
    axis equal
    title(['Final wind a=',string_diag_A])

    figure
    plot_mesh_3d(X,[1,n(1),1,n(2),1,1]), hold on, 
    plot_wind_3d(CX,W,1)
    hold off
    axis equal
    title(['Final wind lowest layer a=',string_diag_A])

    figure
    plot_mesh_3d(X,[1,n(1),1,n(2),1,1]), hold on, 
    plot_wind_3d(CX,W,1:2)
    hold off
    axis equal
    title(['Final wind lowest layers a=',string_diag_A])
    
    figure
    height=1;
    wind_at_h(X,CX,W,[20,20,1],[2,18,1,9,height,height]); hold on
    plot_mesh_3d(X,[1,n(1),1,n(2),1,1])
    hold off
    axis equal
    title(['Final wind with a=',string_diag_A,' at ',num2str(height),' above terrain'])

end

% condition_number=scond(K)
n
