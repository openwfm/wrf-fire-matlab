function path_struct = graph_dets(w,cull)
%w = read_wrfout_tign(wrfout)

[fire_name,save_name,prefix] = fire_choice();
red = subset_domain(w);
time_bounds(2) = red.max_tign;
time_bounds(1) = red.min_tign;

% figures
fig.fig_map=0;
fig.fig_3d=0;
fig.fig_interp=0;
p = sort_rsac_files(prefix);

g_str = 'g.mat';
if ~exist(g_str,'file')
    g = subset_l2_detections(prefix,p,red,time_bounds,fig);
    save(g_str, 'g', '-v7.3');
else
    g = [];
    reload_dets = input_num('Reload detections? 1 = yes',0);
    if reload_dets == 1
        g = subset_l2_detections(prefix,p,red,time_bounds,fig);
        save(g_str, 'g', '-v7.3');
    else
        load g_graph.mat;
    end
end

close all
pts = [];
min_con = 70;
%make unique ignition point 
for i = 1:length(g)% fprintf('Detections collected \n')
    % figure(1),scatter3(pts(:,2),pts(:,1),pts(:,3));
    % title('full scatter')
    if sum(g(i).det(3:5)) > 0
        fires = g(i).conf >= min_con;
        lons = mean(g(i).lon(fires));
        lats = mean(g(i).lat(fires));
        confs = mean(double(g(i).conf(fires)));
        times = g(i).time-0.25;
        pts = [lats,lons,times,confs];
    end
    if norm(pts) > 0
        break
    end
end


for i = 1:length(g)
    
    % don't use times after model end
    if (sum(g(i).det(3:5)) > 0) && (g(i).time < red.end_datenum)
        fires = g(i).conf >= min_con;
        lons = g(i).lon(fires);
        lats = g(i).lat(fires);
        times = g(i).time*ones(size(lons));
        confs = double(g(i).conf(fires));
        pts = [pts;[lats',lons',times',confs']];
    end
end
% fprintf('Detections collected \n')
% figure(1),scatter3(pts(:,2),pts(:,1),pts(:,3));
% title('full scatter')

%prune the data

n_points = pts(1:cull:end,:,:,:);
% figure(2),scatter3(n_points(:,2),n_points(:,1),n_points(:,3)-red.start_datenum);
% title('Partial scatter')

%make edge weights
n = length(n_points);
%adjacency / distance
a = zeros(n,n);
%velocity
v = a;
%%% figure out way to get max_ automatically
max_t = 1.2;

%maybe change later
pts = n_points;

%convert from points in the scattered data of L2 data to nearest 
%neighbors on the fire grid
grid_pts = fixpoints2grid(red,n_points);


E = wgs84Ellipsoid;

%max distance from ignition
max_d = 0;
%error in time of fire from time of detection
time_err = 0.1;
distant_point = 1;
ig_point = [pts(1,1),pts(1,2)];
for i = 1:n
    time = pts(i,3);
    %%% local time for help in figuring out day/night
    locs(i) = local_time(time);
    i_point = [pts(i,1),pts(i,2)];
    %find furthest detection from ignition
    new_d = distance(ig_point,i_point,E);
    if new_d > max_d
        max_d = new_d;
        distant_point = i;
    end
    %distance from all points
    for j = 1:n
        time_diff = max(time_err,pts(j,3)-time);
        %if (time_diff > 0  && time_diff < max_t)
        if time_diff < max_t
            j_point = [pts(j,1),pts(j,2)];
            a(i,j) = distance(i_point,j_point,E);
            %special treatment for detections on same granule
%             if time_diff == time_err
%                 a(i,j) = sqrt(a(i,j));
%             end
            v(i,j) = a(i,j)/time_diff;
        end
    end
end
fprintf('Matrix ready \n');
%scatter ignition and distant point
% figure(2),hold on,zlabel('Days from start'),xlabel('Lon'),ylabel('Lat')
% scatter3(pts(1,2),pts(1,1),pts(1,3)-red.start_datenum,'*r');
% scatter3(pts(distant_point,2),pts(distant_point,1),pts(distant_point,3)-red.start_datenum,'*r');
% hold off

%calculate max velocity of fire
v_std = std(v(:)); %change to v(v>0)?
v_mean = mean(v(:));
%max_v = abs(max_d/(pts(distant_point,3)-pts(1,3)));
max_v = v_mean;
v_mask = v < max_v;
%filter detections too far apart, fix this KLUGE
speed_penalty = 1.2;
fast = 2.0;
slow = 0.4;
speeders = 0;
for i = 1:n
    if (locs(i) > 8) && (locs(i) < 18)
        sp1 = fast;
    else
        sp1 = slow;
    end
    for j = 1:n
        if (locs(j) > 8) && (locs(j) < 18)
            sp2 = fast;
        else
            sp2 = slow;
        end
        speed_penalty = (sp1+sp2)/2;
        %fprintf('speed penalty = %0.2f \n',speed_penalty)
        if (a(i,j) > 0) && (v(i,j) > speed_penalty*max_v)
            a(i,j) = 0;
            speeders = speeders + 1;
        end
    end
end
fprintf('%d speeders removed \n',speeders)


fg = digraph(a);
figure(3),plot(fg);
path_count = 0;
start_pt = 1;
new_points = pts;
for i = 1:1
    %paths to next 20 detections
    for j = 1:cull:n
        distant_point = j;
        %[p,d] = shortestpath(fg,start_pt,distant_point);
        [p,d] = shortestpath(fg,i,j);
        if ~isempty(p)
            path_count = path_count + 1;
            paths(path_count).p = p;
            paths(path_count).d = d;
            %path confidenc is geometric mean of detections in path
            paths(path_count).c = prod(pts(p,4))^(1/length(p));
            %fprintf('%d points in path \n',length(p))
        end
        figure(2),hold on
        %l2 points
        %plot3(pts(p,2),pts(p,1),pts(p,3)-red.start_datenum,':r');
        %grid points
        plot3(grid_pts(p,4),grid_pts(p,3),pts(p,3)-red.start_datenum,'g');
                for k = 2:length(p)-1
                    scatter3(pts(p(k),2),pts(p(k),1),pts(p(k),3)-red.start_datenum,'*r');
                end
        hold off
%         for k = 1:length(p)-1
%             new_pt = ([pts(p(k),1),pts(p(k),2),pts(p(k),3)]+[pts(p(k+1),1),pts(p(k+1),2),pts(p(k+1),3)])/2;
%             new_points = [new_points;new_pt];
%         end
    end
    path_struct.paths = paths;
    path_struct.graph = fg;
    path_struct.distances = a;
    path_struct.speeds = v;
    path_struct.points = double(pts);
    path_struct.red = red;
    path_struct.new_points = new_points;
    path_struct.grid_pts = grid_pts(:,3:4);
    path_struct.idx = uint8(grid_pts(:,1:2));
end
end

function lt = local_time(time)
    shift = 7;
    lt = 24*((time-shift)-floor(time-shift));
end


