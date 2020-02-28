function plot_mesh_3d(X)
n = size(X{1})-1;
hold on
for i=1:prod(n)
    [xi,yi,zi] = ind2sub(n,i);
    x1 = X{1}(xi,yi,zi);
    x2 = X{1}(xi+1,yi,zi);
    y1 = X{2}(xi,yi,zi);
    y2 = X{2}(xi,yi+1,zi);
    z1 = X{3}(xi,yi,zi);
    z2 = X{3}(xi,yi,zi+1);
    if zi == 1
        plot3([x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],[z1,z1,z1,z1,z1],'color','b')
    end
    plot3([x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],[z2,z2,z2,z2,z2],'color','b')
    plot3([x1,x1],[y1,y1],[z1,z2],'color','b')
    plot3([x1,x1],[y2,y2],[z1,z2],'color','b')
    plot3([x2,x2],[y1,y1],[z1,z2],'color','b')
    plot3([x2,x2],[y2,y2],[z1,z2],'color','b')
end
hold off
view(3)
end