function path_struct = tign_try(w,cpd)
%scatter fake data on TIGN cone, try to recover the shape
%cpd clusters per day of data
%[fire_name,save_name,prefix] = fire_choice();
if ~exist('red.mat','file')
    red = subset_domain(w);
    save red.mat red
else
    load red.mat
end
time_bounds(2) = red.max_tign;
time_bounds(1) = red.min_tign;

% figures
fig.fig_map=0;
fig.fig_3d=0;
fig.fig_interp=0;

[m,n] = size(red.tign);
rm = 0.1;
[ig_x,ig_y] = find(red.tign == min(red.tign(:)));
pts = [red.fxlat(ig_x,ig_y),red.fxlong(ig_x,ig_y),min(red.tign(:)),100,100,1];
idx = [ig_x,ig_y];
for k = 1 : 1
    for i = 1:m
        for j = 1:n
            if red.tign(i,j) < red.end_datenum
                if rand < 0.1
                    %pts = [pts;[lats',lons',times',confs',frps',gran']];
                    pts = [pts;[red.fxlat(i,j),red.fxlong(i,j),red.tign(i,j),100,100,2*round(red.tign(i,j)-red.start_datenum)]];
                    idx = [idx;[i,j]];
                end
            end
        end
    end
end
if exist('pts.mat','file')
    clear pts idx~
    load pts.mat
else
    save pts.mat pts idx
end

figure,mesh(red.fxlong,red.fxlat,red.tign-red.start_datenum),hold on
scatter3(pts(:,2),pts(:,1),pts(:,3)-red.start_datenum,'r*')
[st1,st2]=sort(pts(:,3));
pts(:,1) = pts(st2,1);
pts(:,2) = pts(st2,2);
pts(:,3) = pts(st2,3);

cull = 1;
%cluster the data
t1 = red.min_tign;
t2 = red.max_tign;
clusters = round((t2-t1)*cpd);
%clusters = round(length(pts)/20);
%clusters = 20;
[s_idx,s_c] = kmeans(pts(:,1:2),clusters);
%scatter the clusters with coloring
figure,scatter3(pts(s_idx==1,2),pts(s_idx==1,1),pts(s_idx==1,3)-red.start_datenum);
hold on
for i = 2:clusters
  scatter3(pts(s_idx==i,2),pts(s_idx==i,1),pts(s_idx==i,3)-red.start_datenum);
end
hold off

%make edge weights
n = length(pts);
%adjacency / distance matrix
a = zeros(n,n);
%velocity matrix
v = a;
%time matrix
t = a;
%cone volume matrix
%cv = a;
grid_pts = pts(:,1:2);

%distance matrices..
%% computing distance between points using GPS coords
E = wgs84Ellipsoid;

%make cluster center distance matrix
clust_dist=zeros(clusters);
%needed ????
for i = 1:clusters
    i_clust = [s_c(i,1),s_c(i,2)];
    for j = 1:clusters
        j_clust = [s_c(j,1),s_c(j,2)];
        clust_dist(i,j) =  distance(i_clust,j_clust,E);
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
raw_dist = a;

%start filtering distance
cluster_mult = 0.5;
for i = 1:n
    for j = 1:n
        %         % make points in same cluster close
        if (a(i,j) > 0) && (s_idx(i) == s_idx(j)) %in same cluster
            a(i,j) = cluster_mult*a(i,j);
        end
        % make points in different clusters further apart
%         if (a(i,j) > 0) && (s_idx(i)~=s_idx(j))
%             a(i,j) = a(i,j) + clust_dist(s_idx(i),s_idx(j));
%         end
    end
end


%make paths
fg = digraph(a);
%figure(3),plot(fg);
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
        plot3(pts(p,2),pts(p,1),pts(p,3)-red.start_datenum,'g');
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
    
end

path_struct.raw_dist = raw_dist+raw_dist';
path_struct.paths = paths;
path_struct.graph = fg;
path_struct.distances = a;
path_struct.speeds = v;
path_struct.points = double(pts);
path_struct.red = red;
%    path_struct.new_points = new_points;
path_struct.grid_pts = grid_pts;
path_struct.idx = idx;




end


function lt = local_time(time)
    shift = 7;
    lt = 24*((time-shift)-floor(time-shift));
end