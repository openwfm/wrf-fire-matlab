function plot_mesh_3d(X,bbox)
if exist('bbox','var')
    for i=1:3
        X{i}=X{i}(bbox(1):bbox(2),bbox(3):bbox(4),bbox(5):bbox(6));
    end
end
[n(1),n(2),n(3)] = size(X{1});
hold on
for j=1:n(2)
    for k=1:n(3)
        plot3(X{1}(:,j,k),X{2}(:,j,k),X{3}(:,j,k),'color','b'),
    end
end
for i=1:n(1)
    for k=1:n(3)
        plot3(X{1}(i,:,k),X{2}(i,:,k),X{3}(i,:,k),'color','b'),
    end
end
for i=1:n(1)
    for j=1:n(2)
        plot3(squeeze(X{1}(i,j,:)),squeeze(X{2}(i,j,:)),squeeze(X{3}(i,j,:)),'color','b'),
    end
end
hold off
view(3)
xlabel('x'), ylabel('y'), zlabel('z'),
end