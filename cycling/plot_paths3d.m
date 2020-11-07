function plot_paths3d(path_struct,fig_num)
%input path_sruct = graph_dets(w,cull)
%      fig_num - figure to draw on       
%fg = path_struct.graph;
pts = path_struct.points(:,1:3);
pts(:,1:2) = path_struct.grid_pts;
pts(:,3) = (pts(:,3)-floor(min(pts(:,3))));
paths = path_struct.paths;
figure(fig_num)
hold on

for j = 1:length(paths)%length(pts)
    %c = sqrt(path_struct.paths(j).c/100);
    %plot3(pts(paths(j).p,2),pts(paths(j).p,1),pts(paths(j).p,3),'Color',[c c c]);
    plot3(pts(paths(j).p,2),pts(paths(j).p,1),pts(paths(j).p,3));
end
scatter3(pts(:,2),pts(:,1),pts(:,3),'*r')
hold off
grid on,xlabel('Lon'),ylabel('Lat'),zlabel('Time [days since simulation start]')
title('Shortest Paths and Active Fire Detections')
end