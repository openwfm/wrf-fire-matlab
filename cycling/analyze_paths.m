:q
function a = analyze_paths(ps,new_tign)
%looks at the paths in a path struct and figures out some basic statistics
% input   ps - struct, from ps = graph_dets(w,1);
%         new_tign - analysis from squish(ps,1)

%vector for speeds between points in the paths new_tign
p_speeds = [];
%vector for speeds between pts in tigndist = ps.raw_dist(p2,p1);
t_speeds = [];
speeds_count = 0;
%figure(117),contour(ps.red.fxlong,ps.red.fxlat,ps.red.tign);
hold on
%loop through the paths, collect ROS info from each path segment
%extrapolate one poin on the end of each path
for i = 1:length(ps.paths)
    p = ps.paths(i).p;
    %fprintf('Path %d, %d points in path  \n',i,length(p))
    %find speed between points in the paths
    pl = length(p);
    for j = 1:pl
        %pause for plotting
        %pause(0.01)
        if j > 1
            speeds_count = speeds_count + 1;
            %point pair in path segment
            p1 = p(j-1);
            p2 = p(j);
            path_points(speeds_count,1:2)=[p1,p2];
            p_speeds(speeds_count) = ps.speeds(p2,p1)/(24*3600);
            p_fuel_1(speeds_count) = ps.red.nfuel_cat(ps.idx(p1,1),ps.idx(p1,2));
            p_fuel_2(speeds_count) = ps.red.nfuel_cat(ps.idx(p2,1),ps.idx(p2,2));
            %plot path on countour
            %plot(ps.points(p,2),ps.points(p,1),'r','LineWidth',0.5);
            % find these locations on the grid
            %tign at these points
%             t_p1 = ps.red.tign(ps.idx(p1,1),ps.idx(p1,2));
%             t_p2 = ps.red.tign(ps.idx(p2,1),ps.idx(p2,2));
            % time at L2 data recording
            t_p1 = ps.points(j-1,3);
            t_p2 = ps.points(j,3);
            % time at the analysis
            %t_p1 = new_tign(ps.idx(p1,1),ps.idx(p1,2));
            %t_p2 = new_tign(ps.idx(p2,1),ps.idx(p2,2));
            dt = abs(t_p2-t_p1);
            dist = ps.raw_dist(p2,p1);
            %reall, this is the p_speeds
            t_speeds(speeds_count)= (dist/dt)/(24*3600);
            %elevation data
            e_p1 = ps.red.fhgt(ps.idx(p1,1),ps.idx(p1,2));
            e_p2 = ps.red.fhgt(ps.idx(p2,1),ps.idx(p2,2));
            slopes(speeds_count) = (e_p2-e_p1)/dist;

        end
        %scatter(ps.points(p(end),2),ps.points(p(end),1),'*r')
    end
end
hold off
a.path_points = path_points;
a.p_speeds = p_speeds';
a.t_speeds = t_speeds';
a.p_fuel_1 = p_fuel_1';
a.p_fuel_2 = p_fuel_2';
a.matrix(:,1)=p_speeds';
a.matrix(:,2)=p_fuel_1';
a.matrix(:,3)=p_fuel_2';
a.slopes = slopes';

%%%%%% wtf happened %%%%
%% (figure(117),contour(ps.red.fxlong,ps.red.fxlat,ps.red.tign);;

%look at averages
%number by which to multiply standard deviation
std_mult = 0.5;

unq = unique([a.matrix,a.path_points],'rows');
common_fuel1 = mode(unq(:,2));
common_fuel2 = mode(unq(:,3));
fuel_mask = unq(:,2:3) == common_fuel2;
avg_speed_common_fuel = mean([p_speeds(fuel_mask(:,1)),p_speeds(fuel_mask(:,2))]);
figure,histogram([p_speeds(fuel_mask(:,1)),p_speeds(fuel_mask(:,2))]),title('Histogram ROS in common fuel')
%std deviation of ROS in the common fuel type
std_speed_common_fuel = std([p_speeds(fuel_mask(:,1)),p_speeds(fuel_mask(:,2))]);
fast_common_fuel = avg_speed_common_fuel+std_mult*std_speed_common_fuel;
fast_mask = unq(:,1)>fast_common_fuel;
%locations where the most common fuel is and ROS is one standard deviation
%obove mean
fast_fuel_mask = logical([fast_mask.*fuel_mask(:,1),fast_mask.*fuel_mask(:,2)]);
%plot the location of the fast fuel
%figure(21),contourf(ps.red.fxlong,ps.red.fxlat,ps.red.tign,'--k');
figure(21),mesh(ps.red.fxlong,ps.red.fxlat,new_tign);
hold on
title('Locations of high ROS in common fuel')
fast_count = 0;
point_speeds = unq(:,1);
for i =1:length(unq)
    if  unq(i,1) > fast_common_fuel
        fast_count = fast_count+1;
        
        %scatter(ps.points(unq(i,5),2),ps.points(unq(i,5),1),'*r')
        idx_1 = unq(i,4);
        idx_2 = unq(i,5);
        scatter3(ps.points(unq(i,4),2),ps.points(unq(i,4),1),new_tign(ps.idx(idx_1,1),ps.idx(idx_1,2)),'*r')
        scatter3(ps.points(unq(i,5),2),ps.points(unq(i,5),1),new_tign(ps.idx(idx_2,1),ps.idx(idx_2,2)),'*r')
    end
end


end
