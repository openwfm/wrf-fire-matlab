function params=params_defaults
    params.run_fortran=0;
    params.run_matlab=1;
    params.femwind_fortran_test=0;
    params.test_fortran=0;
    params.graphics=2;  % 1=basic, 2=all
    params.expand=1.2;  % exponential grid expansion in the vertical
    params.mesh_top=1000; % if given, ignore params_expand 
    params.sc_all=[1,2]; % mesh refinements for tests at multiple scales 
    params.sc2_all=[1,2,4];  % additional factors for horizonal mesh extent 
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
    params.maxit_coarse=8; % 2 smoothing, coarse, 2 smoothing, coarse, 2 smoothing
    params.coarsest_iter=100; % 0 = direct solver n>0 number of iterations
    params.nsmooth=3; % smoothing iterations before correcton
    params.coarsening='2 linear';
    % params.coarse_P='variational';  
    params.coarse_K='assembly';
    params.P_by_x=1;  % prolongation by geometrically linear interpolation
    params.smoothing='vertical sweeps';
    % params.smoothing='3D red-black';
    params.restol=1e-6;
    params.exact=0; % compare with exact solution to compute error
    params.slice=0.5; % vertical y slice of error to display, 0 to 1
    params.err_slice_fig=12; % figure number for residual slice
    params.res_slice_fig=13; % figure number for error slice
    params.iterations_fig=14; % figure number for iterations progress
    params.maxaspect=3;  % do not coarsen vertically if vertical layer is too thick
    params.minaspect=1/3; % do not coarsen horizontally if the layer is too thin
    params.levels=15;
    params.apply_coarse_boundary_conditions=1;
    params.nsmooth_coarse=2;
    params.save_files=0; % save progress levels=3, workspace=2 params only=1
    params.save_file_prefix='femwind';  
    %Define Streamline Starting Points: Defined in terms of scale*nelem
    params.in_height_stream = [175]; 
    params.time_stream  = 0;
    params.st_contour = 1; %Produce contour planes in streamlines plot 0 off, 1 on
    params.st_quiver = 1; %Produce vectors along streamlines 0 off, 1 on
end
