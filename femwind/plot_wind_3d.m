function plot_wind_3d(X,W,text,level,scale)

all_levels = false;
if ~exist('level','var')
    all_levels = true;
end
if ~exist('scale','var')
    scale = 1;
end
if ~exist('text','var')
    text="wind"
end

if all_levels
    quiver3(X{1}(:),X{2}(:),X{3}(:),W{1}(:),W{2}(:),vz,scale,'LineWidth',2)
    xlabel('x'), ylabel('y'), zlabel('z')
    title(title)
    hold off
else
    n = size(X{1});
    x = X{1}(:,:,level); y = X{2}(:,:,level); z = X{3}(:,:,level);
    vx = reshape(vx,n); vy = reshape(vy,n); vz = reshape(vz,n);
    u = vx(:,:,level); v = vy(:,:,level); w = vz(:,:,level);
    quiver3(x(:),y(:),z(:),u(:),v(:),w(:),scale,'LineWidth',2)
    xlabel('x'), ylabel('y'), zlabel('z')
    title(text)
    hold off
end
end