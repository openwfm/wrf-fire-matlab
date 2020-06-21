function graph_dets(g)

close all
pts = [];
min_con = 8;
for i = 1:length(g)
    if randn > 0 && sum(g(i).det(3:5)) > 0
        fires = g(i).data >= min_con;
        lons = g(i).xlon(fires);
        lats = g(i).xlat(fires);
        times = g(i).time*ones(length(lons),1);
        pts = [pts;[lats,lons,times]];
    end
end
fprintf('Detections collected \n')
figure(1),scatter3(pts(:,2),pts(:,1),pts(:,3));
title('full scatter')

%prune the data
cull = 4;
n_points = pts(1:cull:end,:,:);
figure(2),scatter3(n_points(:,2),n_points(:,1),n_points(:,3));
title('patial scatter')

%make edge weights
n = length(n_points);
a = zeros(n,n);
max_t = 3.5;

%maybe change later
pts = n_points;
E = wgs84Ellipsoid;

%max distance from ignition
max_d = 0;
distant_point = 1;
ig_point = [pts(1,1),pts(1,2)];
for i = 1:n
    time = pts(i,3);
    i_point = [pts(i,1),pts(i,2)];
    new_d = distance(ig_point,i_point,E);
    if new_d > max_d
        max_d = new_d;
        distant_point = i;
    end
    for j = 1:n
        if (pts(j,3)-time) > 0 && (pts(j,3) - time) < max_t
            j_point = [pts(j,1),pts(j,2)];
            a(i,j) = distance(i_point,j_point,E);
        end
    end
end
fprintf('Matrix ready \n');
%scatter ignition and distant point
figure(2),hold on
scatter3(pts(1,2),pts(1,1),pts(1,3),'*r');
scatter3(pts(distant_point,2),pts(distant_point,1),pts(distant_point,3),'*r');
hold off

fg = digraph(a);
figure(3),plot(fg);


[p,d] = shortestpath(fg,1,distant_point);

figure(2),hold on
plot3(pts(p,2),pts(p,1),pts(p,3),'r');
for k = 2:length(p)-1
%     tail = pts(p(k-1),:,:);
%     head = pts(p(k),:,:);
%     arr = tail-head;
    scatter3(pts(p(k),2),pts(p(k),1),pts(p(k),3),'*r');
    %quiver3(tail(2),tail(1),tail(3),arr(2),arr(1),arr(3));
end
hold off
    

end