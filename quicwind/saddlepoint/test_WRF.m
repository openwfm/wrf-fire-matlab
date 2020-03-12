%ifile = 'wrfinput_d01';
%si = nc2struct(ifile,{'XLONG','XLAT','HGT','U','V','W','PH','PHB'},{'DX','DY'});
%ofile = 'wrfout_d01_2010-07-01_01:00:00';
%so = nc2struct(ofile,{'U','V','W'},{});
s = load('wind.mat');

% transformation matrix from wind at x,y,z coordinates to two constant
% winds at the faces
T=[1,0,0;
   1,0,0;
   0,1,0;
   0,1,0;
   0,0,1;
   0,0,1];
% intial wind at the faces 
v0x = (s.u0(1:end-1,:,:)+s.u0(2:end,:,:))/2;
v0y = (s.v0(:,1:end-1,:)+s.v0(:,2:end,:))/2;
v0z = (s.w0(:,:,1:end-1)+s.w0(:,:,2:end))/2;
v0 = T*[v0x(:),v0y(:),v0z(:)]';
v0 = v0(:);
% final wind at the center from WRF
uWx = (s.u(1:end-1,:,:)+s.u(2:end,:,:))/2;
uWy = (s.v(:,1:end-1,:)+s.v(:,2:end,:))/2;
uWz = (s.w(:,:,1:end-1)+s.w(:,:,2:end))/2;
uW = T*[uWx(:),uWy(:),uWz(:)]';
uW = uW(:);
% geopotential height
gh = (s.ph+s.phb)/9.81;

% dimensions and spacings
n = size(v0x);
h = [s.dx,s.dy];

% create meshes from input file
xm = repmat(s.xlong,[1,1,n(3)]);
ym = repmat(s.xlat,[1,1,n(3)]);
zm = (gh(:,:,1:end-1)+gh(:,:,2:end))/2;
Xm = {xm,ym,zm};
[xx,yy,~] = ndgrid(h(1)*(0:n(1)),h(2)*(0:n(2)),0:n(3));
ghm = (gh(1:end-1,1:end-1,:)+gh(2:end,1:end-1,:)+gh(1:end-1,2:end,:)+gh(2:end,2:end,:))/4;
ghm = [ghm(1,:,:);ghm;ghm(end,:,:)];
ghm = [ghm(:,1,:),ghm,ghm(:,end,:)];
zz = ghm;
X = {xx,yy,zz};

% assembly sparse matrices
[A,D,E,B,C,~] = sparse_assembly(X,h);

% solve 
saddle_sparse

% transform resulting fluxes into cartesian winds at the centers
u=E*B*v;

% transform WRF winds at the face to the center 
uWc = E*uW;
v0c = E*v0;

% plots
figure(1)
mesh(xx(:,:,1),yy(:,:,1),zz(:,:,1)) 
hold on
plot_wind_3d(Xm,v0c,'Initial WRF wind',1);
figure(2), 
mesh(xx(:,:,1),yy(:,:,1),zz(:,:,1)) 
hold on
plot_wind_3d(Xm,uWc,'Final WRF wind',1);
figure(3), 
mesh(xx(:,:,1),yy(:,:,1),zz(:,:,1)) 
hold on
plot_wind_3d(Xm,u,'Final mass-consistent wind',1);
figure(4), 
mesh(xx(:,:,1),yy(:,:,1),zz(:,:,1)) 
hold on
plot_wind_3d(Xm,u-uWc,'Final mass-consistent wind',1);

% differences
ux = reshape(u(1:3:end),n); uy = reshape(u(2:3:end),n); uz = reshape(u(3:3:end),n);
figure(5)
mesh(Xm{1}(:,:,1),Xm{2}(:,:,1),ux(:,:,1)-uWx(:,:,1))
xlabel('x'), ylabel('y'), zlabel('z')
title('Differences in x direction')
figure(6)
mesh(Xm{1}(:,:,1),Xm{2}(:,:,1),uy(:,:,1)-uWy(:,:,1))
xlabel('x'), ylabel('y'), zlabel('z')
title('Differences in y direction')
figure(7)
mesh(Xm{1}(:,:,1),Xm{2}(:,:,1),uz(:,:,1)-uWz(:,:,1))
xlabel('x'), ylabel('y'), zlabel('z')
title('Differences in z direction')

max_diff = max(max(max(u-uWc))),
