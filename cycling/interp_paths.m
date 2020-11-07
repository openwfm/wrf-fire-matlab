function new_paths = interp_paths(ps)
%function adds new points along paths
new_paths = ps;

%save easy named variables from struct
pts = ps.grid_pts;
pts(:,3) = ps.points(:,3);
idx = ps.idx;

n = length(ps.paths)
%number of new points on segment
np = 3;
new_points = [];
news = [];
new_idx = ps.idx;
new_grid_pts = ps.grid_pts;
new_points = ps.points;
new_points2 = ps.points;
new_grids = [];

for i = 1:n
    p = ps.paths(i).p;
    new_p = p;
    pl = length(p);
    if pl > 1
        fprintf('\n Making new points, path %d\n',i)
%         if i == 280
%             pause
%         end
        
        for j = 2:pl
            
            %two points in space, p and q
            p_i = idx(p(j-1),1);%+rm*round(randn);
            p_j = idx(p(j-1),2);%+rm*round(randn);
            q_i = idx(p(j),1);%+rm*round(randn);
            q_j = idx(p(j),2);%+rm*round(randn);
            %distance between points in path
            dist = ps.raw_dist(p(j-1),p(j));
            np = round(dist/600)+1;
            fprintf('Segment %d distance: %f new points: %d \n',j-1,dist,np);
            %linear interpolation of new points along the line
            new_lats = linspace(pts(p(j-1),1),pts(p(j),1),np)' ...
                + 1/1000*randn(np,1);
            new_lons = linspace(pts(p(j-1),2),pts(p(j),2),np)' ...
                + 1/1000*randn(np,1);
            new_time = linspace(pts(p(j-1),3),pts(p(j),3),np)';
            %average confidence of each endpoint on segment
            new_conf = mean(ps.points(j-1,4),ps.points(j,4))*ones(np,1);
            new_frps = linspace(ps.points(j-1,4),ps.points(j,4),np)';
            new_gran = ps.points(j-1,6)*ones(np,1);
            news = [news;[new_lats,new_lons,new_time]];
            %new_grids = [new_grids ;fixpoints2grid(ps.red,[new_lats,new_lons])];
            %figure(33),hold on,scatter3(new_points(:,3),new_points(:,4),new_time)
            %new_p = 
            new_points=[new_points;[new_lats(2:end-1),new_lons(2:end-1), ...
                                    new_time(2:end-1),new_conf(2:end-1),new_frps(2:end-1),new_gran(2:end-1)]];
            %new_grid_pts = [new_grid_pts; [new_grids(2:end-1,3),new_grids(2:end-1,4)]];
            %new_idx = [new_idx; [new_grids(2:end-1,1),new_grids(2:end-1,2)]];
            %adding row position of new path points in case its needed
            %later
            idx_end = length(new_points);
            idx_start = idx_end-np+3;
            %new_p = [new_p(1:j-1) idx_start:idx_end new_p(j:end)];
            new_p = [new_p(1:find(new_p==p(j-1))) idx_start:idx_end new_p(find(new_p==p(j)):end)];
            
        end

    end 
    %interpolate the whole path with spline
    new_paths.paths(i).p = new_p;
    if pl > 1
        pplat = csaps(ps.points(p,3),ps.points(p,1),0.95,new_points(new_p,3));
        pplon = csaps(ps.points(p,3),ps.points(p,2),0.95,new_points(new_p,3));
        new_points(new_p,1)=pplat;
        new_points(new_p,2)=pplon;
    end
end % for i
new_grids = fixpoints2grid(ps.red,[new_points(:,1),new_points(:,2)]);
new_paths.grid_pts = [new_grids(:,3),new_grids(:,4)];
new_paths.idx = [new_grids(:,1),new_grids(:,2)];
new_paths.points = new_points;
end % function