function new_path_struct = extrapolate_path_pt(ps,new_tign,w)
%put one point on the end of each path
% input ps - path struct from from ps = graph_dets(w,1);
%       new_tign - analysis from squish(ps,1)
%       w - w = read_wrfout_tign(f), wrfout struct with fxlong,fxlat, tign
%       , etc
%only extrapolate on paths ending within a day of ps.red.end_datenum
new_pt_ct = 0;
%when, beyond red.end_datenum to extraploate to
time_extend = 0.1;
stored_ps = ps;
%vector for locations where fuel need to be dried
fast_list = []
for i = 1:length(ps.paths)
    p = ps.paths(i).p;
    end_pt = p(end);
    if (length(p) > 1) && (ps.red.end_datenum - ps.points(end_pt,3) < 0.25)
        %fprintf('Path length and end time OK \n')
        pt_1 = p(end-1);
        %vector from pt_1 to end_pt
        %calc 3  line equation : r = r_0 + t*v
        v = ps.points(end_pt,1:3)-ps.points(pt_1,1:3);
        %check to make sure time difference is OK
        if v(1,3) > 0
            num_pts = length(ps.points);
            new_pt = num_pts+1;
            new_pt_ct = new_pt_ct +1;
            t_f = ps.red.end_datenum + time_extend;
            r_0 = ps.points(end_pt,1:3);
            t = (t_f-r_0(1,3))/v(1,3);
            %new point r
            r = r_0 + t*v;
            ps.paths(i).p(end+1) = new_pt;
            %new point in list, sam confidence al last point on the path being extended
            ps.points(new_pt,:) = ps.points(end_pt,:);
            ps.points(new_pt,1:3) = r;
            ps.new_points(new_pt,:)=ps.points(new_pt,:);
            %fix new point on the grid
            np = [ps.points(new_pt,1),ps.points(new_pt,2)];
            try
                [new_i,new_j,new_lat,new_lon] = fixpt(ps.red,np);
                ps.idx(new_pt,1)=uint8(new_i);ps.idx(new_pt,2)=uint8(new_j);
                ps.grid_pts(new_pt,1)=new_lon;ps.grid_pts(new_pt,2)=new_lat;
                if ps.points(new_pt,3)-new_tign(new_i,new_j)<1
                    fast_list = [fast_list;new_pt];
                end
                
            catch
                % or just duplicate last point in path
                fprintf('weird point outside of domain or something ..\n')
                ps.idx(new_pt,:)=ps.idx(end_pt,:);
                ps.grid_pts(new_pt,:) = ps.grid_pts(end_pt,:);
                ps.points(new_pt,:) = ps.points(end_pt,:);
                ps.new_points(new_pt,:)=ps.points(end_pt,:);
            end
        end
    end
end
%kluster the data bakfire, head fire
cluster = kmeans(ps.grid_pts(fast_list,:),2);
c1 = fast_list(cluster == 1);
c2 = fast_list(cluster == 2);


fprintf('%d new points were added \n',new_pt_ct)
% figure(1),mesh(ps.red.fxlong,ps.red.fxlat,new_tign)
% time_mask = ps.points(:,3)>ps.red.end_datenum;
% figure(1),hold on,scatter3(ps.points(time_mask,2),ps.points(time_mask,1),ps.points(time_mask,3))
% loop to sort new points where we can change FMC
figure,mesh(ps.red.fxlong,ps.red.fxlat,new_tign)
hold on,scatter3(ps.points(c1,2),ps.points(c1,1),ps.points(c1,3),'b*');
hold on,scatter3(ps.points(c2,2),ps.points(c2,1),ps.points(c2,3),'r*');

ps.fast_list = fast_list;
ps.cluster = cluster;

new_path_struct = ps;

end
