function [dx,dy,nan_msk] = fire_gradients(lon,lat,tign,unit)
%computes the gradients of the fire arrival time cone
% inputs : tign - fire arrival time matrix, lon = fxlon, lat = fxlat
% output:  ux,uy   - vector components of gradients in the x and y directions
% unit = 1 ==> return unit vectors
E = wgs84Ellipsoid;
% [n,m] = size(lon);
% hx = distance(lat(1,round(m/2)),lon(1,1),lat(1,round(m/2)),lon(end,end),E)/m;
% hy = distance(lat(1,1),lon(round(n/2),1),lat(end,end),lon(round(n/2),1),E)/n;


[aspect,slope,dy,dx] = gradientm(lat,lon,tign,E);
%set all gradients at top of cone to NaN
t_top = tign > max(tign(:))-0.1;
dx(t_top) = NaN;
dy(t_top) = NaN;

%[dx,dy] = gradient(tign);
%using Sobel edge detector
% [dx,dy] = imgradientxy(tign);
% dx = dx/8;
% dy = dy/8;
%make unit vectors
if unit == 1
x = dx./sqrt(dx.^2+dy.^2);
y = dy./sqrt(dx.^2+dy.^2);
dx = x;
dy = y;
end

% figure,contour(fxlong,fxlat,tign)
% hold on
% quiver(fxlong,fxlat,ux,uy)
% hold off

%[gy,gx] = imgradientxy(tign);
mx = isnan(dx);
my = isnan(dy);
nan_msk = logical(mx.*my);
%plot the vectors
fprintf('there were %d NaN in the gradient \n',sum(nan_msk(:)));
% figure,quiver(lon(~nan_msk),lat(~nan_msk),dx(~nan_msk),dy(~nan_msk))
% title('Gradient')
% figure,scatter(lon(nan_msk),lat(nan_msk));
% title(['Location of ',num2str(sum(nan_msk(:))),' NaN in gradient'])
end % function