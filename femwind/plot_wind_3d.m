function plot_wind_3d(CX,W,level,scale)
% plot_wind_3d(CX,W,level,scale)
% in:
%     CX     {x,y,z} 3D coordinates where the wind vectors are located
%     W      {u,v,w} 3D wind vectors
%     level  vector of number of vertical levels to display
%            default if none or empty: all levels
%     scale  scaling of wind vector arrow, see doc quiver

all_levels = false;
[n(1),n(2),n(3)]=size(CX{1});
if ~exist('level','var') || isempty(level)
    level=1:n(3);
end
if ~exist('scale','var')
    scale = 0.9;
end

quiver3(CX{1}(:,:,level),CX{2}(:,:,level),CX{3}(:,:,level),...
   W{1}(:,:,level), W{2}(:,:,level), W{3}(:,:,level),...
   scale,'LineWidth',2)
xlabel('x'), ylabel('y'), zlabel('z')
end