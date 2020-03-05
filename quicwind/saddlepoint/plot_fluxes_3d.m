function plot_fluxes_3d(X,U,h,scales)
if ~exist('h','var')
    h = [1,1,1];
end
if ~exist('scales','var')
    max_u = max([U(1:6:end);U(2:6:end)]);
    max_v = max([U(3:6:end);U(4:6:end)]);
    max_w = max([U(5:6:end);U(6:6:end)]);
    max_wind = [max_u,max_v,max_w];
    scales = h./(3*max_wind);
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
    z111 = X{3}(xi,yi,zi);
    z112 = X{3}(xi,yi,zi+1);
    z211 = X{3}(xi+1,yi,zi);
    z212 = X{3}(xi+1,yi,zi+1);
    z121 = X{3}(xi,yi+1,zi);
    z122 = X{3}(xi,yi+1,zi+1);
    z221 = X{3}(xi+1,yi+1,zi);
    z222 = X{3}(xi+1,yi+1,zi+1);
    zu1 = (z111+z112+z121+z122)/4;
    zu2 = (z211+z212+z221+z222)/4;
    zv1 = (z111+z112+z211+z212)/4;
    zv2 = (z121+z122+z221+z222)/4;
    zw1 = (z111+z211+z121+z221)/4;
    zw2 = (z112+z212+z122+z222)/4;
    pp{i} = quiver3([x1,x2,xm,xm,xm,xm],[ym,ym,y1,y2,ym,ym],[zu1,zu2,zv1,zv2,zw1,zw2],[U(s(1:2))',0,0,0,0],[0,0,U(s(3:4))',0,0],[0,0,0,0,U(s(5:6))'],0,'LineWidth',2); 
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