function [ growth_struct] = growth_score( wrfout,prefix )
%function [ gs ] = growth_score( wrfout,prefix )
% inputs:
%   wrfout - string, parth to wrfout file for evaluation
%   prefix - string, path to directory with stellite dat
% outputs:
%   growth_struct- struct with the following:
%      gs - mean error in area gowth rate of forecast vs. area
%       from detections
%      data_time - vector with times of sat observations, (days from sim
%       start)
%      sat_area, fore_area - vectors with estimated areas of fires (grid nodes
%       as units)

% make reduced structure
% if wrfout(end) ~= 't'
%     w = read_wrfout_tign(wrfout);
% else
%     fprintf('wrfout is a mat file, loading \n');
%     load(wrfout);
% end

w = read_wrfout_tign(wrfout);
red = subset_domain(w);
time_bounds = [red.start_datenum red.end_datenum];

% figures
fig.fig_map=0;
fig.fig_3d=0;
fig.fig_interp=0;

%find detections
p = sort_rsac_files(prefix);
if exist('g_full.mat','file')
    load g_full.mat
else
    g = load_subset_detections(prefix,p,red,time_bounds,fig);
    save g_full.mat g;
end


%granules with active fire detections
data_count = 0;

min_confidence = 7;
for i=1:length(g)
    fprintf('Loading granule %d \n',i)
    %
    % add to detections already on fire grid
    
    if i == 1
        fire_mask = g(i).data >=min_confidence;
        sat_lons = g(i).xlon(fire_mask);
        sat_lats = g(i).xlat(fire_mask);
    end
    fire_mask = g(i).data >= min_confidence;
    lon_update = g(i).xlon(fire_mask);
    lat_update = g(i).xlat(fire_mask);
    if sum(fire_mask(:)) > 0
        fprintf('Active fire in domain\n')
        data_count = data_count+1;
        sat_lons = double([sat_lons(:);lon_update(:)]);
        sat_lats = double([sat_lats(:);lat_update(:)]);
        
        fprintf('Computing Satellite area \n')
        %draw polygon around the fire grid detections
        sat_boundary = boundary(sat_lons,sat_lats);
        [sat_in,sat_on] = inpolygon(red.fxlong(:),red.fxlat(:),sat_lons,sat_lats);
        %figure,scatter(
        sat_area(data_count) = sum(sat_in) + sum(sat_on);
        data_time(data_count) = g(i).time;
        
        fprintf('Computing forecast area \n')
        fore_mask = red.tign < g(i).time;
        fore_area(data_count) = sum(fore_mask(:));
        %diff_fore_area(i) = fore_area(i) - fore_area(i-1)
        
        
        fprintf('Computing Growth Rates \n')
        if data_count == 1
            diff_sat_area(data_count) = sat_area(data_count);
            diff_fore_area(data_count) = fore_area(data_count);
        else
            diff_sat_area(data_count) = sat_area(data_count)-sat_area(data_count-1);
            diff_fore_area(data_count) = fore_area(data_count)-fore_area(data_count-1);
        end
        %else diff_sat_area(1) = temp_sat_area
        
        fprintf('Computing relative error \n')
        %also plot the sequences of growth rates in same axes   
        end
end

%convert data_time to days since simulation
data_time = (data_time - red.start_datenum);
if strcmp(prefix,'../TIFs/')
    term = 40;
else
    term = length(data_time);
end

fprintf('%d granules had active fire data\n',data_count);
figure,plot(data_time(1:term),sat_area(1:term))
hold on, plot(data_time(1:term),fore_area(1:term))
legend('Area within detection perimeter','Area within forecast perimeter')
xlabel('Simulation Time [days]')
ylabel('Simulation Area [grid nodes]')
title(wrfout);
hold off

% figure,plot(data_time(1:term),diff_sat_area(1:term))
% hold on, plot(data_time(1:term),diff_fore_area(1:term))
% legend('Change in IR Perimter Area','Change in Forecast Area')
% xlabel('Simulation Time [days]')
% ylabel('Change in Area [grid nodes]')
% hold off
gs = mean(abs(sat_area-fore_area));

%compute rates
sat_rate = sat_area;
fore_rate = fore_area;
for i = 2:length(sat_area)
    sat_rate(i) = sat_area(i) - sat_area(i-1);
    fore_rate(i) = fore_area(i) - fore_area(i-1);
end


growth_struct.gs = gs;
growth_struct.data_time = data_time;
growth_struct.sat_area = sat_area;
growth_struct.fore_area = fore_area;
growth_struct.sat_rate = sat_rate;
growth_struct.fore_rate = fore_rate;
end

