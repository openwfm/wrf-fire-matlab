function [dx,dy] = fire_gradients(fxlong,fxlat,tign,unit)
%computes the gradients of the fire arrival time cone
% inputs : tign - fire arrival time matrix
% output:  ux,uy   - unit vector components of gradients in the x and y directions
% unit = 1 ==> return unit vectors
lon = fxlong;
lat = fxlat;
E = wgs84Ellipsoid;
% [n,m] = size(lon);
% hx = distance(lat(1,round(m/2)),lon(1,1),lat(1,round(m/2)),lon(end,end),E)/m;
% hy = distance(lat(1,1),lon(round(n/2),1),lat(end,end),lon(round(n/2),1),E)/n;


[aspect,slope,dy,dx] = gradientm(lat,lon,tign,E);
%[dx,dy] = gradient(tign);
%using Sobel edge detector
% [dx,dy] = imgradientxy(tign);
% dx = dx/8;
% dy = dy/8;
%make unit vectors
if unit == 1
dx = dx./sqrt(dx.^2+dy.^2);
dy = dy./sqrt(dx.^2+dy.^2);
end

% figure,contour(fxlong,fxlat,tign)
% hold on
% quiver(fxlong,fxlat,ux,uy)
% hold off

%[gy,gx] = imgradientxy(tign);


end % function