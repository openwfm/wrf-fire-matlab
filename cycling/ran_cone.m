function rc = ran_cone(w)
%makes cone structure
close all
f = @(x,y,x0,y0,m) m*sqrt((x-x0).^2+(y-y0).^2);

[p,q] = size(w.fxlong);
tx = linspace(-10,10,p);
ty = linspace(-10,10,q);
[x,y] = meshgrid(tx,ty);

x0 = 0;
y0 = 0;
z0 = 0;
z_max = 20;
radius = 9;

r = sqrt(x0^2+y0^2);
% m <= (z_max-z0)/(9-r)
m = (1/2+rand)*(z_max-z0)/(radius-r);
%m = 1;

z = f(x,y,x0,y0,m);
z(z>z_max) = z_max;

n = round(2*z_max);

base = 2*rand;
for i = 1:n
    %vertext of new cone
    z0 = i/n;
    z0 = 1/2*z_max*rand;
%     new_x0 = new_x0+9/n*randn;
%     new_y0 = new_y0+9/n*randn;
    new_x0 = 2.5*i/n*randn;
    new_y0 = 2.5*i/n*randn;
    
    r = sqrt(new_x0^2+new_y0^2);
    new_m = (z_max-z0)/(radius+2*randn-r)+base+randn;
    c(i,1) = new_x0;
    c(i,2) = new_y0;
    z2 = f(x,y,new_x0,new_y0,new_m);

%     noise = 10/z_max*randn(p,q);
%     z2 = z2+noise;
    
    alpha = 1/2;
    z2(z2>z_max) = z_max;
    %figure,mesh(x,y,z),hold on,mesh(x,y,z2),hold off;
    %z(z2>z0) = max(z(z2>z0),z2(z2>z0));
    z(z2>z) = alpha*z(z2>z)+(1-alpha)*z2(z2>z);
    z = imgaussfilt(z,1/n);
    %z = smooth_up(x,y,z);
    z(z>z_max) = z_max;
    
    %figure(234),contour(x,y,z,20,'k')
    
end  %for loop    

%z = smooth_up(x,y,z);
%figure,mesh(x,y,z);
%quick_mesh(z)
%figure,plot(c(:,1),c(:,2));
figure,contour(x,y,z,20,'k');


rc = z;

end % function
