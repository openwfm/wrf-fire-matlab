function [wz] = perim_cone(n,w)

[x,y,z] = perim();  %in-linefunction below
for i = 2:n
    [x1,y1,z1] = perim();
    z = min(z,z1);
end
close all
%figure,contourf(x,y,z,20,'k');
%figure,mesh(x,y,z);

%time dimensions
t_min = min(w.tign_g(:));
t_max = max(w.tign_g(:));


%compute fire area
min_lat = min(w.fxlat(:));max_lat = max(w.fxlat(:));
min_lon = min(w.fxlong(:));max_lon = max(w.fxlong(:));
E = wgs84Ellipsoid;
dlon= distance(min_lat,min_lon,min_lat,max_lon,E);
dlat= distance(min_lat,min_lon,max_lat,min_lon,E);

area_mask = z<max(z(:));
fire_area = dlon*dlat*sum(area_mask(:))/1e5;
fprintf('Firea area: %f ha \n',fire_area)
area_day = fire_area/(t_max-t_min)*3600*24; %m^2/day
fprintf('Burn rate: %f   ha/day\n',area_day);



%fit cone shape to wrfout size
[p,q] = size(w.fxlong);
%[m,n] = size(x);
x_min = min(x(:));
x_max = max(x(:));
y_min = min(y(:));
y_max = max(y(:));
tx = linspace(x_min,x_max,q);
ty = linspace(y_min,y_max,p);
[qx,qy] = meshgrid(tx,ty);
qz = interp2(x,y,z,qx,qy);


z_min = min(qz(:));
z_max = max(qz(:));

wz = t_min + (t_max-t_min)/(z_max-z_min)*(qz-z_min);
t_msk = wz<max(wz(:));

%plot
figure,contourf(w.fxlong,w.fxlat,wz,20,'k');
% minlon = min(w.fxlong(t_msk));
% maxlon = max(w.fxlong(t_msk));
% dx = maxlon - minlon;
% minlat = min(w.fxlat(t_msk));
% maxlat = max(w.fxlat(t_msk));
% dy = maxlat - minlat;
% xlim([minlon-dx/4 maxlon+dx/4])
% ylim([minlat - dy/4 maxlat + dy/4])

end


function [x,y,z] = perim()

fprintf('Making Cone\n')
close all
l = 24;
theta = linspace(0,2*pi,l);
% ellipse paramters
a = 2+5*rand;
b = 2+5*rand;
%regular coordinates of ellipse
x = a*cos(theta);
y = b*sin(theta);

%figure(23),scatter(x,y)
for i = 1:l
    x(i) = x(i)+randn/i;
    y(i) = x(i)+randn/i;
    for j = 1:l
       idx = mod(j,l)+1;
       x(idx) = 1/2*(x(j)+x(idx))+2*randn/a;
       y(idx) = 1/2*(y(j)+y(idx))+2*randn/b;
    end
%     figure(23),scatter(x,y)
%     pause(0.1)
end
% pause(1)
%close all
%recenter
xbar = mean(x);
ybar = mean(y);
x=x-xbar;
y=y-ybar;
% figure(23),scatter(x,y)
% hold on
% for i = 1:l
%     plot([0 x(i)],[0 y(i)])
%     pause(0.2)
% end
% hold off

t_final = 10;
pts = [0 0 0];
for i = 1:l
    xpts = linspace(0,x(i),t_final)+1/8*randn(1,t_final);
    ypts = linspace(0,y(i),t_final)+1/8*randn(1,t_final);
    zpts = linspace(0,t_final,t_final)+1/2*randn(1,t_final);
    pts = [pts;[xpts',ypts',zpts']];
end

% figure,scatter3(pts(:,1),pts(:,2),pts(:,3))
dim1 = 75;
dim2 = 75;
F = scatteredInterpolant(pts(:,1),pts(:,2),pts(:,3),'linear','nearest');
d = round(max(max(pts(:,1:2))))+15+10*rand;
[x,y] = meshgrid(linspace(-d,d,dim1),linspace(-d,d,dim2));
z = F(x,y);

for i = 1:2
%alpha = rand;
z1 = smooth_up(x,y,z);
z2 = imgaussfilt(z,3+rand);
%z = alpha*z1+(1-alpha)*z2;
z = 1/2*(z1+z2);
z(z>t_final) = t_final;
end

msk = logical(ones(size(z)));
bdry = 10;
msk(bdry:end-bdry,bdry:end-bdry)=false;
tmax = min(z(msk));
z(z>tmax) = tmax;

%figure,mesh(x,y,z)
%figure,contour(x,y,z,20,'k');


end
