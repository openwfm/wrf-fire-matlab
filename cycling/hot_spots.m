function sp = hot_spots(p_gran,red)
%function lookts at satellie data between perimeter observations and
%clusters data into back fire , fire froint
%inputs:
%p_gran - struct with perimeters in the format of l2 active fire granules
%red = red = subset_domain(w)


%display and chose which perimeters to work with
pl = length(p_gran);
for i = 1:pl
    fprintf('%d :%s  %s \n',i,p_gran(i).file,datestr(p_gran(i).time))
end
pt = input_num('Which two perimeters to use as time points? Enter in time order.',[1 2])
t1 = p_gran(pt(1)).time;
t2 = p_gran(pt(2)).time;
fprintf('Start date: %s\n',datestr(t1))
fprintf('End date: %s\n',datestr(t2))
%subset detection data
[~,~,prefix,perim] = fire_choice();
time_cushion = 0.4;
time_bounds(1) = t1-time_cushion;
time_bounds(2) = t2+time_cushion;
fig.fig_map=0;
fig.fig_3d=0;
fig.fig_interp=0;
p = sort_rsac_files(prefix);
%load satellite data
g_str = 'g_hot.mat';
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
        load g_hot.mat;
    end
end
%allow for ealier and later data
min_con = 70;
pts = [];
for i = 1:length(g)
    if (sum(g(i).det(3:5)) >0) && (g(i).time > t1-time_cushion) && (g(i).time < t2+time_cushion)
        fprintf('reading satellite data')
        fires = g(i).conf >= min_con;
        lons = g(i).lon(fires);
        lats = g(i).lat(fires);
        times = g(i).time*ones(size(lons));
        confs = double(g(i).conf(fires));
        frps = g(i).power(fires);
        gran = i*ones(size(lons));
        pts = [pts;[lats',lons',times',confs',frps',gran']];
    else
        fprintf('No Fires, or fires outside of time boundary. Time = %s\n',datestr(g(i).time))
    end
end
%plot perims and the detection data
figure,scatter(p_gran(pt(1)).lon,p_gran(pt(1)).lat,'k')
hold on
scatter(p_gran(pt(2)).lon,p_gran(pt(2)).lat,'b')
scatter(lons,lats,'*r')
legend('Perimeter 1','Perimeter 2','Fire Detections')
end %function