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
red = subset_domain(w,1);


%use the wrfout file to find subset of detections

time_bounds(2) = red.max_tign;
%time_bounds(2) = 7.354637292824074e+05
time_bounds(1) = red.min_tign;
%time_bounds(1) = time_bounds(2)-3;
fire_choice = input_num('which fire? Patch: [0], Camp: [1], Cougar: [3]',1)
cycle = input_num('Which cycle? ',0)
if fire_choice == 1
    prefix='../campTIFs/';
elseif fire_choice == 0
    prefix='../TIFs/';
else
    prefix = '../cougarTIFs/'
end
det_list=sort_rsac_files(prefix);
fig.fig_map=0;
fig.fig_3d=0;
fig.fig_interp=0;
if exist('g_match.mat','file')
    load('g_match.mat')
else
    g = load_subset_detections(prefix,det_list,red,time_bounds,fig);
    save g_match.mat g
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
    if max(lat_update(:)) > 46.2 | min(lat_update(:)) < 46.05
        fprintf('bad granule : %s \n',g(i).file)
    end
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
decimate = 1;
fire_lon = fire_lon(1:decimate:end);
fire_lat = fire_lat(1:decimate:end);
%scatter(fire_lon,fire_lat,'*','r')
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
title({'Satellite Fire Detections and Forecast Perimeter',score_str});
legend({'Forecast perimeter','Detections inside perimeter','Detections outside perimeter'});
%ylim([39.5 39.95])
%xlim([-121.9 -121.3])
hold off

% plots simulation perimeter
% scatter(fire_lon(fire_on),fire_lat(fire_on),'r*')


    



end

