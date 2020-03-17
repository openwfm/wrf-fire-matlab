function [ plotters ] = plot_full_sim(fig_num, fire, wrfout, title_string )
% inputs:
%   fig_num - integer, figure number
%   fire - string, either "Patch" or "Camp"
%   wrfout - path to a wrfout file
%   title_string, string, title of figure
%plots fire cone

% make reduced structure
w = read_wrfout_tign(wrfout);
red = subset_domain(w);

% get detection data
if fire(1) == 'P'
    prefix = '../TIFs/';
else 
    prefix = '../campTIFs/';
end

% figures
fig.fig_map=0;
fig.fig_3d=0;
fig.fig_interp=0;


time_bounds = [red.start_datenum red.end_datenum];
p = sort_rsac_files(prefix);
if ~exist('g_full.mat','file')
    g = load_subset_detections(prefix,p,red,time_bounds,fig);
    save g_full.mat g
else
    load g_full.mat
end

if fire(1) == 'P'
    plot_state(fig_num,red,title_string,red.tign,g,time_bounds(1:2))
else
    red.fxlong = red.fxlong(1:10:end,1:10:end);
    red.fxlat = red.fxlat(1:10:end,1:10:end);
    plot_state(fig_num,red,title_string,red.tign(1:10:end,1:10:end),g,time_bounds(1:2))
    figure(fig_num)
    zlim([0 2.2]);
    %zlim([min(red.tign(:))/(24*3600) max(red.tign(:))/(24*3600)])
end

end

