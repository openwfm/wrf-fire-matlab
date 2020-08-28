function tign = squish2(ps)
%function takes path structure and makes new estimated tign
%inputs  ps = graph_dets(w) or ps = cluster_paths(w), struct with points, paths, etc.

%copy variables to shorter names
pts = ps.grid_pts;
pts(:,3) = ps.points(:,3);
new_pts = pts;
tign = ps.red.tign;
idx = ps.idx;
new_idx = idx;
fignum = 111;

%make new points by interpolation along the paths 
for i = 1:length(ps.paths)
    p = ps.paths(i).p;
    for j = 1:length(p)
%         p_i = idx(p(j),1);
%         p_j = idx(p(j),2);
        % new point
        if j > 1
            w1 = 0.5;
            w2 = 0.5;
%             frp1 = ps.points(p(j-1),5);
%             frp2 = ps.points(p(j),5);
%             w1 = frp1/(frp1+frp2);
%             w2 = frp2/(frp1+frp2);
            new_lat = w1*pts(p(j-1),1)+w2*pts(p(j),1);
            new_lon = w1*pts(p(j-1),2)+w2*pts(p(j),2);
            new_t = w1*pts(p(j-1),3)+w2*pts(p(j),3);
            [new_i,new_j,new_lat,new_lon]= fixpt(ps.red,[new_lat,new_lon]);
            new_pts=[new_pts;[new_lat new_lon new_t]];
            new_idx = [new_idx;[new_i,new_j]];
        end
    end % for j  
end %for i
%new points ready to work with
[new_pts,u1,u2]=unique(new_pts,'rows');
% will this uniques command give good results used separately ???
idx = new_idx(u1,1:2);
times = new_pts(:,3);
lats = new_pts(:,1);
lons = new_pts(:,2);

%sort data by time
[st1,st2] = sort(times);
times = times(st2);
lats = lats(st2);
lons=lons(st2);
idx(:,1:2) = idx(st2,1:2);

n = length(new_pts);
end_time = max(max(times),ps.red.end_datenum);
start_time = min(times);
total_time = ceil(end_time-start_time);

tign = end_time*ones(size(ps.red.fxlong));
temp_tign = tign;
pts_tign = tign;
lons_set = [];
lats_set = [];
times_set = [];

%make random patch around detections
rm = 0;
for i = 1:n;% total_time %half day interval
    if rm ~= 0
        pts_tign(idx(i,1)-round(rm*rand):idx(i,1)+round(rm*rand),idx(i,2)-round(rm*rand):idx(i,2)+round(rm*rand))=times(i);
    end
    %draw polygon around every 20 detections
    if mod(i,20) == 0 || i == n
    pt_set = times <= times(i);
    lons_set = [lons_set;lons(pt_set)];
    lats_set = [lats_set;lats(pt_set)];
    times_set = [times_set;times(pt_set)];
%     figure,scatter3(lons_set,lats_set,times_set)
    %convex hull operations
%     dt = delaunayTriangulation(lats_set,lons_set);
%     c = convexHull(dt);
%     [in,on] = inpolygon(ps.red.fxlat,ps.red.fxlong,lats_set(c),lons_set(c));
%     figure,plot(lons_set(c),lats_set(c));
    [in,on] = inpolygon(ps.red.fxlat,ps.red.fxlong,lats_set,lons_set);
    %in = logical(in+on);
    temp_tign(in) = max(times_set);
    temp_tign(~in) = end_time;
    tign = min(tign,temp_tign);
    %tign = imgaussfilt(tign,1);
    figure(fignum)
    mesh(ps.red.fxlong,ps.red.fxlat,tign)
    end
end
%tign=min(tign,pts_tign);
tign = imgaussfilt(tign,2/3,'FilterDomain','frequency');
fprintf('Blurring ... \n')
figure(fignum),mesh(ps.red.fxlong,ps.red.fxlat,tign)
xlabel('Lon'),ylabel('Lat'),zlabel('Time'),title('Interpolated by Polygons')
%resize, if small matrices were used
if isfield(ps.red,'red')
    if size(ps.red.fxlong) ~= size(ps.red.red.fxlong)
        fprintf('resizing to original size of grid \n');
        [n,m] = size(ps.red.red.tign);
        la = linspace(ps.red.red.min_lat,ps.red.red.max_lat,m);
        lo = linspace(ps.red.red.min_lon,ps.red.red.max_lon,n);
        [lat,lon] = meshgrid(la,lo);
        tign = griddata(ps.red.fxlat(:),ps.red.fxlong(:),tign(:),lat,lon);
    end
end


end
