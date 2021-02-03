function params=femwind_test(params)
% usage:
% params=femwind_test    output params structure and exit
% femwind_test(params)   run with params specified and copy to output
% 
disp('femwind_test')
% test components of the femwind system

if ~exist('params','var') | isempty(params)
    params.graphics=1;  % 1=basic, 2=all
    params.expand=1.2;  % exponential grid expansion in the vertical
    params.mesh_top=1000; % if given, ignore params_expand 
    params.sc_all=[1,2]; % mesh refinements for tests at multiple scales 
    params.sc2_all=[1,2,4,8,16,32];  % additional factors for horizonal mesh extent 
    params.nelem3=[50,50,10]; % base size in elements, horizontal=2*odd 
    params.h=[30,30,2]; % base mesh spacing before scaling
    params.a=[1 1 1]; % penalty factors in x y z directions
    params.initial_wind='log'; % or uniform
    params.roughness_height=0.5;
    params.terrain_shape='hill'; % terrain for add_terrain_to_mesh
    params.terrain_top='shift'; % mesh top treatment for add_terrain_to_mesh
    params.terrain_height=0.1; % terrain height as part of domain height
    params.solver='2-level' ; % see sparse_solve.m
    params.maxit=50; % max iterations
    params.coarsening='2 linear';
    params.smoothing='vertical sweeps';
    % params.smoothing='3D red-black';
    params.nsmooth=3; % smoothing iterations before correcton
    params.restol=1e-6;
    params.exact=0; % compare with exact solution to compute error
    params.slice=0.5; % vertical y slice of error to display, 0 to 1
    params.err_slice_fig=12; % figure number for residual slice
    params.res_slice_fig=13; % figure number for error slice
    params.iterations_fig=14; % figure number for iterations progress
    params.maxaspect=3;  % do not coarsen vertically if vertical layer is too thick
    params.minaspect=1/3; % do not coarsen horizontally if the layer is too thin
    params.levels=3;
    params.apply_coarse_boundary_conditions=1;
    params.nsmooth_coarse=2;
    params.maxit_coarse=8; % 2 smoothing, coarse, 2 smoothing, coarse, 2 smoothing
    params.save_files=2; % save progress
    return 
end

params

if isfield(params,'mesh_top')
    if params.mesh_top>0
        disp('given params.mesh_top>0, computing params.expand') 
        nz = params.nelem3(3); % elements in the vertical direction
        a = params.mesh_top/params.h(3); % desired height as multiple of first layer
        params.expand = findq(a,nz);
    end
end
params

for sc = params.sc_all
    for sc2 = params.sc2_all
        nel = sc*params.nelem3;  % elements in the 3 directions
        nel(1:2)=nel(1:2)*sc2
        h = params.h/sc;
        fprintf('mesh of %ix%ix%i cells\n',nel(1),nel(2),nel(3))
        params.id=sprintf('%ix%ix%i',nel); % to pass around 
        string_diag_A=sprintf('%g %g %g',params.a); % for figure titles
        A = diag(params.a.^2);
        lambda = zeros(prod(nel+1),1); % placeholder solution

        % creating the grid
        expand=params.expand;
        X = regular_mesh(nel,h,params.expand^(1/sc));
        X = add_terrain_to_mesh(X,...
            params.terrain_shape,params.terrain_top,params.terrain_height);
        [CX,CH] = center_mesh(X); % get midpoints of elements

        % initial wind at the centers of the elements
        rng(1);
        switch params.initial_wind
            case 'uniform'
                disp('initial wind uniform in x direction')
                U0={ones(nel),zeros(nel),zeros(nel)};
            case 'random-z'
                % to test iterative methods with non-smooth initial error
                disp('initial wind uniform in x direction random in z direction')
                U0={ones(nel),zeros(nel),randn(nel)};
            case 'random-xz'
                % to test iterative methods with non-smooth initial error
                disp('initial wind uniform in x direction and z direction')
                U0={ones(nel),zeros(nel),ones(nel)};
            case 'log'
                disp('initial wind log profile in x direction')
                U0={log(max(1,CH/params.roughness_height)),zeros(nel),zeros(nel)};
            otherwise
                error(['unknown initial wind ',params.initial_wind])
        end
        if params.graphics>0
            disp('graphics: problem setup')
            % show mesh
            figure(1),clf
            plot_mesh_3d(X)
            axis equal
            title('The wind mesh, wind vector in centers, lambda in corners')
        end

        if params.graphics>0
            % show initial wind
            figure(2),clf
            plot_mesh_3d(X,[1,nel(1)+1,1,nel(2)+1,1,1]), hold on, 
            plot_wind_3d(CX,U0)
            hold off
            axis equal
            title('Initial wind')
        end

        if params.graphics>1
            % show initial wind
            figure(3),clf
            plot_mesh_3d(X,[1,nel(1),1,nel(2)+1,1,1]), hold on, 
            plot_wind_3d(CX,U0,1)
            hold off
            axis equal
            title('Initial wind lowest layer')
        end

        % assemble sparse system matrix
        [K,F,~] = sparse_assembly(A,X,U0,lambda,params);

        % dirichlet boundary conditions
        [K,F]=apply_boundary_conditions(K,F,X);

        % solve the equations
        % [lambda,it] = sparse_solve(K,F,X,'s');
        [lambda,it,rate(sc,sc2),XC,P] = sparse_solve(K,F,X,params);

        % assemble final wind
        [~,~,W] = sparse_assembly(A,X,U0,lambda,params);

        if params.graphics>1
            disp('graphics: solution')

            % plot resulting wind
            figure(4),clf
            plot_mesh_3d(X,[1,nel(1)+1,1,nel(2)+1,1,1]), hold on, 
            plot_wind_3d(CX,W)
            hold off
            axis equal
            title(['Final wind a=',string_diag_A])

            figure(5),clf
            plot_mesh_3d(X,[1,nel(1)+1,1,nel(2)+1,1,1]), hold on, 
            plot_wind_3d(CX,W,1)
            hold off
            axis equal
            title(['Final wind lowest layer a=',string_diag_A])

            figure(6),clf
            plot_mesh_3d(X,[1,nel(1)+1,1,nel(2)+1,1,1]), hold on, 
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
            plot_mesh_3d(X,[1,nel(1)+1,1,nel(2)+1,1,1])
            hold off
            axis equal
            title(['Final wind with a=',string_diag_A,' at ',num2str(height),' above terrain'])
            
            figure(8),clf
            wind_streamlines(X,W,params)
        end    
        if params.save_files>0,
            save -v7.3
        end
    end
end
