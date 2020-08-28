function path_struct = cluster_paths(w,cull)
% assign shortest paths, using clustering
% inputs - w = read_wrfout_tign(f)
%          cull - number fo using smaller data sets
% output - cp , struct with path info

[fire_name,save_name,prefix] = fire_choice();
red = subset_domain(w);
%shrink the size for large matrices
target_size = 200;
if length(red.tign) > target_size

    [m,n] = size(red.tign);
    %shrink_factor
    %sf = 4;
    max_dim = max(m,n);
    sf = target_size/max_dim;
    n =round(n*sf);
    m = round(m*sf);
    red_copy = red;
    red = subset_small(red,m,n);
end
time_bounds(2) = red.max_tign;
time_bounds(1) = red.min_tign;

% figures
fig.fig_map=0;
fig.fig_3d=0;
fig.fig_interp=0;
p = sort_rsac_files(prefix);

%load satellite data
g_str = 'g_cluster.mat';
if ~exist(g_str,'file')
    %loading L2 data
    g = subset_l2_detections(prefix,p,red,time_bounds,fig);
    save(g_str, 'g', '-v7.3');
else
    g = [];
    reload_dets = input_num('Reload detections? 1 = yes',0);
    if reload_dets == 1
        g = subset_l2_detections(prefix,p,red,time_bounds,fig);
        save(g_str, 'g', '-v7.3');
    else
        load g_cluster.mat;
    end
end

close all
pts = [];
%minimum detection confidence level
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
        frps = mean(g(i).power(fires));
        gran = i;
        pts = [lats,lons,times,confs,frps,gran];
    end
    if norm(pts) > 0
        break
    end
end

%can change end time for comparisons
end_time = red.end_datenum;
for i = 1:length(g)
    % don't use times after model end
    if (sum(g(i).det(3:5)) > 0) && (g(i).time < end_time)
        fires = g(i).conf >= min_con;
        lons = g(i).lon(fires);
        lats = g(i).lat(fires);
        times = g(i).time*ones(size(lons));
        confs = double(g(i).conf(fires));
        frps = g(i).power(fires);
        gran = i*ones(size(lons));
        pts = [pts;[lats',lons',times',confs',frps',gran']];
    end
end

%prune the data
n_points = pts(1:cull:end,:,:,:,:,:);

%% for computing distance between points using GPS coords
% also used for finding aspect of the slope, for clustering
E = wgs84Ellipsoid;
[aspect,slope,dy,dx] = gradientm(red.fxlat,red.fxlong,red.fhgt,E);
clst_pts =  fixpoints2grid(red,n_points);
ig_pt = [mean(clst_pts(:,3)),mean(clst_pts(:,4))];
for i = 1:length(clst_pts)
   pt_1 = [ig_pt(1,1),clst_pts(i,4)];
   pt_2 = [clst_pts(i,3),clst_pts(i,4)];
   d_lon = sign(clst_pts(1,4)-clst_pts(i,4))*distance(ig_pt,pt_1,E);
   d_lat = sign(clst_pts(1,3)-clst_pts(i,3))*distance(pt_2,pt_1,E);
   cp(i,:) = [clst_pts(i,:),d_lat,d_lon];
    %% work out x-y coordinate with pt 1 as origin
end


%cluster the data 
dt = 3*ceil(g(end).time - g(1).time);
space_clusters = dt; %days
[s_idx,s_c] = kmeans(pts(:,1:2),space_clusters);
%clustering using aspect, not good
%[s_idx,s_c] = kmeans(cp(:,5:6),space_clusters);

% find optimal cluster k
% max_clusts = 20;
% s_pts = pts(:,1:2);
% klist=2:max_clusts;%the number of clusters you want to try
% myfunc = @(X,K)(kmeans(X, K));
% eva = evalclusters(s_pts,myfunc,'CalinskiHarabasz','klist',klist)
% classes=kmeans(s_pts,eva.OptimalK);
% dt = eva.OptimalK;

% %spatial clusters scatter plot
% figure,scatter(pts(s_idx==1,1),pts(s_idx==1,2));
% hold on% dt = round(g(end).time - g(1).time);% dt = round(g(end).time - g(1).time);
% space_clusters = dt; %days
% [s_idx,s_c] = kmeans(pts(:,1:2),space_clusters);


% for i = 2:dt
%   scatter(pts(s_idx==i,1),pts(s_idx==i,2));
% end
% hold off

%scatter 3d
figure,scatter3(pts(s_idx==1,2),pts(s_idx==1,1),pts(s_idx==1,3));
hold on
for i = 2:dt
  scatter3(pts(s_idx==i,2),pts(s_idx==i,1),pts(s_idx==i,3));
end
hold off

%make edge weights
n = length(n_points);
%adjacency / distance matrix
a = zeros(n,n);
%velocity matrix
v = a;
%time matrix
t = a;
%cone volume matrix
cv = a;
%%% figure out way to get max_t automatically
% maximum allowed time between nodes in the graph to allow them to be
% connected
max_t = 1.9*(24*3600);

%maybe change later
pts = n_points;

%convert from points in the scattered data of L2 data to nearest 
%neighbors on the fire grid
grid_pts = fixpoints2grid(red,n_points);

%% computing distance between points using GPS coords
E = wgs84Ellipsoid;

%make cluster center distance matrix
clust_dist=zeros(space_clusters);
for i = 1:space_clusters
    i_clust = [s_c(i,1),s_c(i,2)];
    for j = 1:space_clusters
        j_clust = [s_c(j,1),s_c(j,2)];
        clust_dist(i,j) =  sqrt((i_clust(1)-j_clust(1))^2+((i_clust(2)-j_clust(2))^2));
    end
end

%max distance from ignition
max_d = 0;
%error in time of fire from time of detection
time_err = 0.2;
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
        time_diff = (pts(j,3)-time)*(24*3600);

        t(i,j) = time_diff;  %pts(j,3)-time;
        if time_diff > 0
            j_point = [pts(j,1),pts(j,2)];
            a(i,j) = distance(i_point,j_point,E);
            v(i,j) = a(i,j)/max(time_diff,0.1);
        end
    end
end

%fix up triangular matrices
t_mask = t >= 0;
% a = a+a';
raw_dist = a;
fprintf('matrices done\n')

%start filtering distance
cluster_mult = 0.5;
for i = 1:n
    for j = 1:n
%         % make points in same cluster close
        if (a(i,j) > 0) && (s_idx(i) == s_idx(j)) %in same cluster
            a(i,j) = cluster_mult*a(i,j);
        end
        % make points in different clusters further apart
        % if (a(i,j) > 0) && (s_idx(i)~=s_idx(j))
        %    a(i,j) = a(i,j) + clust_dist(s_idx(i),s_idx(j));
        % end
    end
end


%make paths
fg = digraph(a);
figure(3),plot(fg);
path_count = 0;

% new_points = pts;
start_pt = 1;
for i = 1:start_pt
    for j = 1:cull:n
        % finds shortest path between points i,j
        % p is the points in the path, d is the total distance
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
        for k = 1:length(p)
            scatter3(pts(p(k),2),pts(p(k),1),pts(p(k),3)-red.start_datenum,'*r');
        end
        hold off
        %         % add a new point to the list by interpolation
        %         for k = 1:length(p)-1
        %             new_pt = ([pts(p(k),1),pts(p(k),2),pts(p(k),3)]+[pts(p(k+1),1),pts(p(k+1),2),pts(p(k+1),3)])/2;
        %             new_points = [new_points;new_pt];
        %         end
    end
    path_struct.raw_dist = raw_dist+raw_dist';
    path_struct.paths = paths;
    path_struct.graph = fg;
    path_struct.distances = a;
    path_struct.speeds = v;
    path_struct.points = double(pts);
    path_struct.red = red;
    %    path_struct.new_points = new_points;
    path_struct.grid_pts = grid_pts(:,3:4);
    path_struct.idx = grid_pts(:,1:2);
end




end %function

function lt = local_time(time)
    shift = 7;
    lt = 24*((time-shift)-floor(time-shift));
end
