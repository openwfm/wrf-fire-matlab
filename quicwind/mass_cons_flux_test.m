function err = mass_cons_flux_test
disp('mass_cons_flux_test - testing mass consistent flux approximation with terrain')

% hh=rand(1,3);
% hh=[1,1,1];
% [x,y,z]   = ndgrid(hh(1)*[0:nx],hh(2)*[0:ny],hh(3)*[0:nz]);
% X={x,y,z};
mesh_len=[50,50,50];
mesh_disp=[20,20,20];
%h=rand(1,3);
h=ones(1,3);
X = regular_mesh(mesh_len,h,1.2);
X = add_terrain_to_mesh(X,'hill','squash',0.05);
x = X{1};
y = X{2};
z = X{3};

U0 = grad3z(ones(size(x)-1),mesh_len);
u0 = 1;
hz = z-z(:,:,1);
h0 = .5;
h_at_u = (hz(:,1:end-1,1:end-1)+hz(:,2:end,1:end-1)+hz(:,1:end-1,2:end)+hz(:,2:end,2:end))/4;
logh = log(max(h_at_u));
U0{1}=u0*log(max(h_at_u/h0,1))/max(logh(:));
U0{2}=0*ones(size(U0{2}));
U0{3}=0*ones(size(U0{3}));
%{ 
% still under development (initial wind following terrain and smoothly go horizontal)
[nx,ny,nz] = size(X{1});
h_at_v = (hz(1:end-1,:,1:end-1)+hz(2:end,:,1:end-1)+hz(1:end-1,:,2:end)+hz(2:end,:,2:end))/4;
hmax = min(min(hz(:,:,end)));
theta = zeros(nx-1,ny,nz); % theta on the ground level at midpoints in the x direction only
theta(:,:,1) = atan((z(2:end,:,1)-z(1:end-1,:,1))./(x(2:end,:,1)-x(1:end-1,:,1)));
for i=2:size(U0{3},3),
    theta(:,:,i) = max(1-h_at_v(:,:,i-1)/hmax,0).*theta(:,:,1);
end
theta_at_u = (theta(1:end-1,:,:)+theta(2:end,:,:))*.5;
U0{1}(2:end-1,:,:) = U0{1}(2:end-1,:,:).*cos(theta_at_u); 
U0{2}=0*ones(size(U0{2}));
U0{3}=U0{3}.*sin(theta);
%}
plot_wind_above_terrain(U0,X,mesh_disp)

[U,Lambda_d] = mass_cons_flux(U0,X,'direct','check'); 
%plot_wind_above_terrain(U,X,mesh_disp)
[V,Lambda_pcg,err] = mass_cons_flux(U0,X,'pcg','check');
plot_wind_above_terrain(V,X,mesh_disp)
direct_vs_pcg = big(cell2vec(U)-cell2vec(V))
lambda_diff = big(Lambda_d-Lambda_pcg)

% plot_wind_above_terrain(U,X,[1 1 1])
end