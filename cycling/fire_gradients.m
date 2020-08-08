function [dx,dy] = fire_gradients(fxlong,fxlat,tign,unit)
%computes the gradients of the fire arrival time cone
% inputs : tign - fire arrival time matrix
% output:  ux,uy   - unit vector components of gradients in the x and y directions
% unit = 1 ==> return unit vectors

[dx,dy] = gradient(tign);
%make unit vectors
if unit == 1
dx = dx./sqrt(dx.^2+dy.^2);
dy = dy./sqrt(dx.^2+dy.^2);
end

% figure,contour(fxlong,fxlat,tign)
% hold on
% quiver(fxlong,fxlat,uy,ux)
% hold off

%[gy,gx] = imgradientxy(tign);


end % function