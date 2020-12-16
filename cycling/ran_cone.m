function rc = ran_cone(w)
%makes cone structure
close all
f = @(x,y,x0,y0,m) m*sqrt((x-x0).^2+(y-y0).^2)+1.5*cos(4*sqrt((x-x0).^2+(y-y0).^2))-1;

[p,q] = size(w.fxlong);
tx = linspace(-10,10,p);
ty = linspace(-10,10,q);
[x,y] = meshgrid(tx,ty);

x0 = randn;
y0 = randn;
z0 = 0;
z_max = 20;
radius = 9+2*rand;

r = sqrt(x0^2+y0^2);
% m <= (z_max-z0)/(9-r)
%m = (1/2+rand)*(z_max-z0)/(radius-r);
m = (1+2*rand)*(z_max-z0)/(radius-r);
%m = 1;

z = f(x,y,x0,y0,m);
z(z>z_max) = z_max;

n = round(1/2*z_max);

base = rand;
new_x0 = x0;
new_y0 = y0;
theta = 2*pi*rand;
for i = 1:n
    %vertext of new cone
%     z0 = i/n;
%     z0 = 1/4*z_max*rand;
    z0 = z0+rand;
    if mod(i,2)
        new_x0 = (new_x0+rand)*cos(theta);
        new_y0 = (new_y0+rand)*sin(theta);
    else
        new_x0 = 2*i/n*randn;
        new_y0 = 2*i/n*randn;
    end
    
    r = sqrt(new_x0^2+new_y0^2);
    new_m =(1+4*rand)*(z_max-z0)/(radius+2*randn-r)+base+rand;
    c(i,1) = new_x0;
    c(i,2) = new_y0;
    z2 = f(x,y,new_x0,new_y0,new_m);

%     noise = 10/z_max*randn(p,q);
%     z2 = z2+noise;
    
    alpha = 3/4;
    z2(z2>z_max) = z_max;
    %figure,mesh(x,y,z),hold on,mesh(x,y,z2),hold off;
    %z(z2>z0) = max(z(z2>z0),z2(z2>z0));
    z(z2>z) = alpha*z(z2>z)+(1-alpha)*z2(z2>z);
    %z = imgaussfilt(z,1/2);
    %z = smooth_up(x,y,z);
    z(z>z_max) = z_max;
    
    %figure(234),contour(x,y,z,20,'k')
    
end  %for loop    

%z = smooth_up(x,y,z);
%figure,mesh(x,y,z);
%quick_mesh(z)
%figure,plot(c(:,1),c(:,2));
%figure,contour(x,y,z,20,'k');


rc = z;

end % function
