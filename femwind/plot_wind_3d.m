function plot_wind_3d(CX,W,level,scale,stride)
% plot_wind_3d(CX,W,level,scale)
% in:
%     CX     {x,y,z} 3D coordinates where the wind vectors are located
%     W      {u,v,w} 3D wind vectors
%     level  vector of number of vertical levels to display
%            default if none or empty: all levels
%     scale  scaling of wind vector arrow, see doc quiver
%     stride display wind vectors

all_levels = false;
[n(1),n(2),n(3)]=size(CX{1});
if ~exist('level','var') || isempty(level)
    level=1:n(3);
end
if ~exist('scale','var') || isempty(scale)
    scale = 0.9;
end

if ~exist('stride','var') || isempty(stride)
    stride = 1;
end

quiver3(CX{1}(1:stride:end,1:stride:end,level),...
    CX{2}(1:stride:end,1:stride:end,level),...
    CX{3}(1:stride:end,1:stride:end,level),...
    W{1}(1:stride:end,1:stride:end,level),...
    W{2}(1:stride:end,1:stride:end,level),...
    W{3}(1:stride:end,1:stride:end,level),...
    scale,'LineWidth',2)
xlabel('x'), ylabel('y'), zlabel('z')
end