disp('femwind_test')
% test components of the femwind system

if ~exist('params','var')
    disp('params do not exist yet, setting')
    params.graphics=2;  % 1=basic, 2=all
    params.expand=1.2;  % exponential grid expansion in the vertical
    params.sc=[1]; % mesh refinements for tests at multiple scales 
    params.sc2=[1,2,4,8,16];  % additional factors for horizonal mesh extent 
    params.nelem3=[20,20,8]; % base size in elements in the 3 directions 
    params.h=[10,10,10]; % base mesh spacing before scaling
    params.da=[1 1 1]; % penalty factors in x y z directions
    params.initial_wind=1;
    params.terrain_shape='hill'; % terrain for add_terrain_to_mesh
    params.terrain_top='squash'; % mesh top treatment for add_terrain_to_mesh
    params.terrain_height=0.2; % terrain height as part of domain height
    params.solver='2-level' ; % see sparse_solve.m
    params.maxit=50; % max iterations
    params.coarsening='2 linear';
    params.smoothing='vertical lines';
    params.nsmooth=5; % smoothing iterations before correcton
    params.restol=1e-6;
    params.exact=1; % compare with exact solution to compute error
    params.slice=0.5; % vertical y slice of error to display, 0 to 1
    params.err_slice_fig=12; % figure number for residual slice
    params.res_slice_fig=13; % figure number for error slice
end
params

for sc = params.sc
    for sc2 = params.sc2
        n = sc*params.nelem3;  % elements in the 3 directions
        n(1:2)=n(1:2)*sc2
        h = params.h/sc;
        fprintf('mesh of %ix%ix%i cells\n',n(1),n(2),n(3))
        string_diag_A=sprintf('%g %g %g',params.da); % for figure titles
        A = diag(params.da);
        lambda = zeros(prod(n+1),1); % placeholder solution

        % creating the grid
        expand=params.expand;
        X = regular_mesh(n,h,params.expand^(1/sc));
        X = add_terrain_to_mesh(X,...
            params.terrain_shape,params.terrain_top,params.terrain_height);
        CX = center_mesh(X); % get midpoints of elements

        % initial wind at the centers of the elements
        rng(1);
        switch params.initial_wind
            case 1
                disp('initial wind uniform in x direction')
                U0={ones(n),zeros(n),zeros(n)};
            case 2
                % to test iterative methods with non-smooth initial error
                disp('initial wind uniform in x direction random in z direction')
                U0={ones(n),zeros(n),randn(n)};
            case 3
                % to test iterative methods with non-smooth initial error
                disp('initial wind uniform in x direction and z direction')
                U0={ones(n),zeros(n),ones(n)};
                
        end
        if params.graphics>0
            disp('graphics: problem setup')
            % show mesh
            figure(1),clf
            plot_mesh_3d(X)
            axis equal
            title('The wind mesh, wind vector in centers, lambda in corners')
        end

        if params.graphics>1
            % show initial wind
            figure(2),clf
            plot_mesh_3d(X,[1,n(1),1,n(2),1,1]), hold on, 
            plot_wind_3d(CX,U0)
            hold off
            axis equal
            title('Initial wind')

            % show initial wind
            figure(3),clf
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
        % [lambda,it] = sparse_solve(K,F,X,'s');
        [lambda,it,rate(sc,sc2),XC] = sparse_solve(K,F,X,params);

        % assemble final wind
        [~,~,W] = sparse_assembly(A,X,U0,lambda);

        if params.graphics>1
            disp('graphics: solution')

            % plot resulting wind
            figure(4),clf
            plot_mesh_3d(X,[1,n(1),1,n(2),1,1]), hold on, 
            plot_wind_3d(CX,W)
            hold off
            axis equal
            title(['Final wind a=',string_diag_A])

            figure(5),clf
            plot_mesh_3d(X,[1,n(1),1,n(2),1,1]), hold on, 
            plot_wind_3d(CX,W,1)
            hold off
            axis equal
            title(['Final wind lowest layer a=',string_diag_A])

            figure(6),clf
            plot_mesh_3d(X,[1,n(1),1,n(2),1,1]), hold on, 
            plot_wind_3d(CX,W,1:2)
            hold off
            axis equal
            title(['Final wind lowest layers a=',string_diag_A])

        end

        if params.graphics>0
            disp('graphics: wind_at_h')

            figure(7),clf
            height=10;
            bbox=[min(X{1}(:)),max(X{1}(:)),...
                min(X{2}(:)),max(X{2}(:)),...
                height,height];
            [XH,WH]=wind_at_h(X,CX,W,[20,20,1],bbox);
            plot_wind_3d(XH,WH)
            hold on
            plot_mesh_3d(X,[1,n(1),1,n(2),1,1])
            hold off
            axis equal
            title(['Final wind with a=',string_diag_A,' at ',num2str(height),' above terrain'])
            
            figure(8),clf
            wind_streamlines(X,W,params)
        end    
    end
end
