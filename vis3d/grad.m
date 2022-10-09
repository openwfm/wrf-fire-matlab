function [dfdx,dfdy]=grad(x,y,f)
% compute directional central derivatives of 2D field
% in:
%    x, y   curvilinear rectangular grid of x,y nodal coordinates 
%    f      values of function on the nodes
% out:
%    dfdx,dfdy   derivatives of f along the x,y lines by central differences 

dxx=  x(3:end, 2:end-1)-x(1:end-2, 2:end-1);
dyx=  y(3:end, 2:end-1)-y(1:end-2, 2:end-1);
dx=sqrt(dxx.^2+dyx.^2);
dfdx=(f(3:end, 2:end-1)-f(1:end-2, 2:end-1))./dx;

dxy=  x(2:end-1, 3:end)-x(2:end-1, 1:end-2);
dyy=  y(2:end-1, 3:end)-y(2:end-1, 1:end-2);
dy=sqrt(dxy.^2+dyy.^2);
dfdy=(f(2:end-1, 3:end)-f(2:end-1, 1:end-2))./dy;

end