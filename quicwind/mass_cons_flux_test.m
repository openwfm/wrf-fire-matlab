function err = mass_cons_flux_test
disp('mass_cons_flux_test - testing mass consistent flux approximation with terrain')

%% Settings
% options:
%   'constant' : constant wind to x direction, 
%   'log' : vertical log profile,
init_wind = 'log'; 
terrain_following = true; % still under developmentssss
u0 = 1; % initial wind magnitude
h0 = .5; % roughness high: only for log profile cases ('log' or 'terrain')
mesh_len=[20,20,20];
mesh_disp=[20,20,20];
h=ones(1,3);
X = regular_mesh(mesh_len,h,1.2);
X = add_terrain_to_mesh(X,'hill','squash',0.05);
x = X{1};
y = X{2};
z = X{3};
[nx,ny,nz] = size(X{1});

%% Initial wind
% declare initial wind
U0 = grad3z(ones(size(x)-1),mesh_len);
switch init_wind
    case 'constant'
        U0{1} = u0*ones(size(U0{1}));
        U0{2} = 0*ones(size(U0{2}));
        U0{3} = 0*ones(size(U0{3}));
    case 'log'   
        hz = z-z(:,:,1);
        h_at_u = (hz(:,1:end-1,1:end-1)+hz(:,2:end,1:end-1)+hz(:,1:end-1,2:end)+hz(:,2:end,2:end))/4;
        mlogh = log(max(h_at_u(:)));
        U0{1} = u0*log(max(h_at_u/h0,1))/mlogh;
        U0{2} = 0*ones(size(U0{2}));
        U0{3} = 0*ones(size(U0{3}));
    otherwise
        error('mass_cons_flux_test: unknown init_wind')
end
plot_wind_above_terrain(U0,X,mesh_disp)
if terrain_following
    U0 = terrain_correction(U0,X);
    plot_wind_above_terrain(U0,X,mesh_disp)
end

%% Direct solve
[U,Lambda_d] = mass_cons_flux(U0,X,'direct','check'); 
plot_wind_above_terrain(U,X,mesh_disp)

%% PCG solve
[V,Lambda_pcg,err] = mass_cons_flux(U0,X,'pcg','check');
plot_wind_above_terrain(V,X,mesh_disp)

%% Test results
direct_vs_pcg = big(cell2vec(U)-cell2vec(V))
lambda_diff = big(Lambda_d-Lambda_pcg)

end