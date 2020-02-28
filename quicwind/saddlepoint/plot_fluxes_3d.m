function plot_fluxes_3d(X,U,h,scales)
if ~exist('h','var')
    h = [1,1,1];
end
if ~exist('scales','var')
    max_u = max([U(1:6:end);U(2:6:end)]);
    max_v = max([U(3:6:end);U(4:6:end)]);
    max_w = max([U(5:6:end);U(6:6:end)]);
    max_wind = [max_u,max_v,max_w];
    scales = h./(2*max_wind);
    scales(max_wind < eps) = eps;
end
plot_mesh_3d(X), 
hold on
n = size(X{1})-1;
factor = 6;
pp = cell(1,prod(n));
for i=1:prod(n)
    s=(i-1)*factor+1:i*factor; 
    [xi,yi,zi]=ind2sub(n,i); 
    x1 = X{1}(xi,yi,zi);
    x2 = X{1}(xi+1,yi,zi);
    xm = (x1+x2)/2;
    y1 = X{2}(xi,yi,zi);
    y2 = X{2}(xi,yi+1,zi);
    ym = (y1+y2)/2;
    z1 = X{3}(xi,yi,zi);
    z2 = X{3}(xi,yi,zi+1);
    zm = (z1+z2)/2;
    pp{i} = quiver3([x1,x2,xm,xm,xm,xm],[ym,ym,y1,y2,ym,ym],[zm,zm,zm,zm,z1,z2],[U(s(1:2))',0,0,0,0],[0,0,U(s(3:4))',0,0],[0,0,0,0,U(s(5:6))'],0,'LineWidth',2); 
end
% scales
su = scales(1); sv = scales(2); sw = scales(3);
for i=1:prod(n)
    hu = get(pp{i},'UData');
    hv = get(pp{i},'VData');
    hw = get(pp{i},'WData');
    set(pp{i},'UData',su*hu,'VData',sv*hv,'WData',sw*hw);
end
hold off,
xlabel('x'), ylabel('y'), zlabel('z')
end