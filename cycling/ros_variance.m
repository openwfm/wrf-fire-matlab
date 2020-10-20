function [ros,ros_mean] = ros_variance(n)
%estimates variance in ros between two detections
%n - number of random points to use in simulation

real_data = input_num('Use real data? 0 = no 1 = yes',1,1);

if real_data == 0
    %coordinates of detections
    d1 = [0 0 0];
    d2 = [3 4 12];
    
    %reported
    dist = sqrt((d2(1)-d1(1))^2+(d2(2)-d1(2))^2);
    time = d2(3) - d1(3);
    ros = dist/time;
    fprintf('ROS = %f\n',ros)
    
    %statistics
    sig1 = 1;
    sig2 = 1;
    time_diff = 0;
    
else
    %real data
    load p_struct.mat
    ps = p;
    clear p
    num_paths = length(ps.paths);
    path_num = max(2,round(num_paths*rand));
    path = ps.paths(path_num);
    path_length = length(path);
    path_end = max(2,round(path_length*rand));
    p1 = path.p(path_end -1);
    p2 = path.p(path_end);
    fprintf('Generating stats for path %d, points %d and %d \n',path_num,path_end-1,path_end);
    lon1 = ps.points(p1,2);
    lat1 = ps.points(p1,1);
    time1 = ps.points(p1,3);
    lon2 = ps.points(p2,2);
    lat2 = ps.points(p2,1);
    time2 = ps.points(p2,3);
    
    E = wgs84Ellipsoid;
    dist = distance(lat1,lon1,lat2,lon2,E);
    time = time2-time1;
    ros = dist/(time*24*3600);
    fprintf('ROS = %f\n',ros)
    
    %coordinates of detections
    d1 = [0 0 0];
    d2 = [0 dist time];
    
    %statistics
    sig1 = 1000;
    sig2 = 1000;
    time_diff = 0.25;
end %if real_data == 0



%n = 200;
r1 = sig1*randn(n,3)+d1;
r1(:,3) = d1(3)-time_diff*rand(n,1);
r2 = sig2*randn(n,3)+d2;
r2(:,3) = d2(3)-time_diff*rand(n,1);

dr = r2-r1;
dist_rand = sqrt(dr(:,1).^2+dr(:,2).^2);
time_rand = dr(:,3);
if real_data == 0
    ros_rand = dist_rand./time_rand;
else
    ros_rand = dist_rand./(24*3600*time_rand);
end
if n < 1000
    figure,scatter3(r1(:,1),r1(:,2),r1(:,3),'*b')
    xlabel('x'),ylabel('y'),zlabel('time')
    title('Random paths between two detections')
    hold on
    scatter3(r2(:,1),r2(:,2),r2(:,3),'*r')
    for i = 1:n
        plot3([r1(i,1) r2(i,1)],[r1(i,2) r2(i,2)],[r1(i,3) r2(i,3)],':')
    end
    hold off
end

figure,histogram(ros_rand)
if real_data == 0
    title('Histogram of ROS')
else
    t_str = sprintf('Histogram of ROS \n Path: %d Points: %d and %d',path_num,path_end-1,path_end);
    title(t_str)
xlabel('ROS'),ylabel('Number')
ros_mean = mean(ros_rand);
fprintf('Mean random ROS: %f\n',mean(ros_rand))


end