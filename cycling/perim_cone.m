function [wz] = perim_cone(n,w)

[x,y,z] = perim(); %in-line function below
for i = 2:n
    [x1,y1,z1] = perim();
    z = min(z,z1);
end
close all
figure,contour(x,y,z,20,'k');
%figure,mesh(x,y,z);

%fit cone shape to wrfout size
[p,q] = size(w.fxlong);
[m,n] = size(x);
x_min = min(x(:));
x_max = max(x(:));
y_min = min(y(:));
y_max = max(y(:));
tx = linspace(x_min,x_max,q);
ty = linspace(y_min,y_max,p);
[qx,qy] = meshgrid(tx,ty);
qz = interp2(x,y,z,qx,qy);

t_min = min(w.tign_g(:));
t_max = max(w.tign_g(:));
z_min = min(qz(:));
z_max = max(qz(:));

wz = t_min + (t_max-t_min)/(z_max-z_min)*(qz-z_min);

end


function [x,y,z] = perim()

fprintf('Making Cone\n')
close all
l = 24;
theta = linspace(0,2*pi,l);
% ellipse paramters
a = 5+10*rand;
b = 5+10*rand;
%regular coordinates of ellipse
x = a*cos(theta);
y = b*sin(theta);

%figure(23),scatter(x,y)
for i = 1:l
    x(i) = x(i)+randn/i;
    y(i) = x(i)+randn/i;
    for j = 1:l
       idx = mod(j,l)+1;
       x(idx) = 1/2*(x(j)+x(idx))+4*randn/a;
       y(idx) = 1/2*(y(j)+y(idx))+4*randn/b;
    end
    %figure(23),scatter(x,y)
end
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
    xpts = linspace(0,x(i),t_final)+1/5*randn(1,t_final);
    ypts = linspace(0,y(i),t_final)+1/5*randn(1,t_final);
    zpts = linspace(0,t_final,t_final)+1/5*randn(1,t_final);
    pts = [pts;[xpts',ypts',zpts']];
end

% figure,scatter3(pts(:,1),pts(:,2),pts(:,3))
dim1 = 75;
dim2 = 70;
F = scatteredInterpolant(pts(:,1),pts(:,2),pts(:,3),'linear','nearest');
d = round(max(max(pts(:,1:2))))+6;
[x,y] = meshgrid(linspace(-d,d,dim1),linspace(-d,d,dim2));
z = F(x,y);

for i = 1:2
z1 = smooth_up(x,y,z);
z2 = imgaussfilt(z,1+rand);
z = 1/2*(z1+z2);
%z(z>t_final) = t_final;
end

msk = logical(ones(size(z)));
bdry = 10;
msk(bdry:end-bdry,bdry:end-bdry)=false;
tmax = min(z(msk));
z(z>tmax) = tmax;

%figure,mesh(x,y,z)
%figure,contour(x,y,z,20,'k');


end
