function [ score ] = match_detections( wrfout, wrf_time_step )
%function score = match_detections( input_args )
%Function evaluates a simulation by finding how well it was able to predict
%where satellites detections would indicate a fire was burning
%Inputs:
%  wrfout: string with path to wrfout file containing the fire arrival
%      time variable tign
%  wrf_time-step: optional string with time step in wrfout to be read
%Output:
%  score: evaluation of the goodness of the fit

%read the wrfout file
% use which time step?
if nargin > 2
    w = read_wrfout_tign(wrfout,wrf_time_step);
else
    if wrfout(end) ~= 't'
        w = read_wrfout_tign(wrfout);
    else
        fprintf('wrfout is already a .mat file, loading \n')
        load(wrfout)
    end
end
red = subset_domain(w);


%use the wrfout file to find subset of detections

time_bounds(2) = red.max_tign;
end_date = datestr(time_bounds(2));
%time_bounds(2) = 7.354637292824074e+05
time_bounds(1) = red.min_tign;
%time_bounds(1) = time_bounds(2)-3;
fire_choice = input_num('which fire? Patch: [0], Camp: [1], Cougar: [3]',1)
cycle = input_num('Which cycle? ',0)
if fire_choice == 1
    fire_name = 'Camp fire';
    save_name = 'camp';
    prefix='../campTIFs/';
elseif fire_choice == 0
    fire_name = 'Patch Springs fire';
    save_name = 'patch';
    prefix='../TIFs/';
else
    fire_name = 'Cougar Creek fire';
    save_name = 'cougar';
    prefix = '../cougarTIFs/';
end

fig.fig_map=0;
fig.fig_3d=0;
fig.fig_interp=0;
g_str = sprintf('g_match_%d.mat',cycle);
if exist(g_str,'file')
    load(g_str)
else
    det_list=sort_rsac_files(prefix);
    g = load_subset_detections(prefix,det_list,red,time_bounds,fig);
    save(g_str, 'g', '-v7.3')
end
%find list of detections
min_con = 7;
for i = 1:length(g)
    if i == 1
        fire_mask = g(i).data >=min_con;
        lons = g(i).xlon(fire_mask);
        lats = g(i).xlat(fire_mask);
    end
    fire_mask = g(i).data >= min_con;
    lon_update = g(i).xlon(fire_mask);
    lat_update = g(i).xlat(fire_mask);
    %looking for bad cougar granule
    %if max(lat_update(:)) > 46.2 | min(lat_update(:)) < 46.05
    %    fprintf('bad granule : %s \n',g(i).file)
    %end
    if sum(fire_mask(:)) > 0
        lons = [lons(:);lon_update(:)];
        lats = [lats(:);lat_update(:)];
    end
    %scatter(lons,lats);
    
end %for i...

%find polygon containing fire perimter from wrf
fire_area = red.tign <= time_bounds(2)-0.001;
% cells within simulation perimeter
fire_lon = red.fxlong(fire_area);
fire_lat = red.fxlat(fire_area);
fire_boundary = boundary(fire_lon,fire_lat);
x = fire_lon(fire_boundary);
y = fire_lat(fire_boundary);


%[fire_in, fire_on] = inpolygon(red.fxlong(:),red.fxlat(:),x,y);
[in, on] = inpolygon(lons,lats,x,y);
%perim_lon = fire_lon(fire_on);
%perim_y = fire_lat(fire_on);
%calculate score


score = (sum(in)+sum(on))/numel(lats);
score_str = sprintf('Cycle %d: %f percent of detections within perimeter',cycle,score*100);
figure
hold on
plot(x,y,'r')
scatter(lons(in),lats(in),'g*')
scatter(lons(~in),lats(~in),'b*')
% t1 = 
% t2 =
% t3 =
title_str = sprintf('Satellite Fire Detections and Forecast Perimeters \n %s \n %s, %s',score_str,fire_name,end_date);
title(title_str)
legend({'Forecast perimeter','Detections inside perimeter','Detections outside perimeter'});
%ylim([39.5 39.95])
%xlim([-121.9 -121.3])
save_str = sprintf('match_%s_%s_%s.fig',save_name,num2str(time_bounds(2)),num2str(cycle))
fprintf('Saving figure as: %s \n',save_str)
savefig(save_str);
saveas(gcf,[save_str '.png']);
hold off

end

