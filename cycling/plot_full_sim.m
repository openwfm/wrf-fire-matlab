function [ plotters ] = plot_full_sim(wrfout )
% inputs:
%   fig_num - integer, figure number
%   fire - string, either "Patch" or "Camp"
%   wrfout - path to a wrfout file
%   title_string, string, title of figure
%plots fire cone

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

fig_num = input_num('Figure number? ',113,1);
% make reduced structure
w = read_wrfout_tign(wrfout);
red = subset_domain(w);


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

title_string = sprintf('%s , cycle %i',fire_name,cycle);
plot_state(fig_num,red,title_string,red.tign,g,time_bounds(1:2))
save_str = sprintf('%s_full_sim_cyc_%d',save_name,cycle)
savefig(save_str);
saveas(gcf,[save_str '.png']);

end

