function X=regular_mesh(n,h,v)
% X=regular_mesh(n,h)
% Creates a regular mesh that is stretched in the vertical direction but
% uniform horizationally.
% Use plot_mesh_3d(X) to see the result.
% in:
%     n           size 3 number of cells in each direction
%     h           size 3 step size in each direction
%     v           size 1 vertical stretch factor
% out:
%     X           cell array with x y z grid coordinates
% 
if ~exist('v','var')
    v=1;
end
zz = zeros(1,n(3)+1);
for i=1:n(3)
    zz(i+1) = zz(i) + h(3) * v^i;
end
[x,y,z] = ndgrid(h(1)*[0:n(1)],h(2)*[0:n(2)],zz);
X = {x,y,z};
fprintf('regular_mesh nodes %i %i %i step %g %g %g expansion %g size %g %g %g\n',...
    n,h,v,n(1)*h(1),n(2)*h(2),zz(end))