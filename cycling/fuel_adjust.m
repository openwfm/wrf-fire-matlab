function fuel_adjust(restart_file)
%function adjust fuel moisture in a  wrf file used for input wrfinput or
%wrfrst

fire_choice = input_num('which fire? Patch: [0], Camp: [1], Cougar: [3]',1)
cycle = input_num('Which cycle? ',0)
load_str = sprintf('p_%d.mat',cycle)
load(load_str);

%nned to beef up what's in red, below
w  = read_wrfout_tign(restart_file);
%use analysis to et region for changing fuel moisture
w.tign_g = p.analysis;
red = subset_domain(w);

%use the wrfout file to find subset of detections

time_bounds(2) = red.max_tign;
time_bounds(1) = red.min_tign;

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
if exist('g_match.mat','file')
    load('g_match.mat')
else
    det_list=sort_rsac_files(prefix);
    g = load_subset_detections(prefix,det_list,red,time_bounds,fig);
    save('g_match.mat', 'g', '-v7.3')
end
%find list of detections
min_con = 7;
for i = 1:length(g)
    if i == 1
        fm = g(i).fxdata >=min_con;
    end
    fm = fm + (g(i).fxdata >=min_con); 
end %for i...

% flatten fm to height 1
fm = fm > 0;

forecast_area = red.tign <= g(i).time;

mult = input_num('Enter fuel moisture multiplier: ',1)
%mult=1.4;
too_far = mult*(forecast_area-fm);
new2far = ones(size(too_far));
too_far_mask = too_far > 0;
new2far(too_far_mask) = too_far(too_far_mask);

fast_fuel = red.nfuel_cat(too_far_mask);
%figure,histogram(fast_fuel(:))

%go to w from red
w_mask = ones(size(w.tign_g));
w_mask(red.ispan,red.jspan) = new2far;


fprintf('%s %s\n','Will write modified time into     ',restart_file)
rewrite_bak=[restart_file,'.bak'];
q=input_num(['1 to copy ',restart_file,' to ',rewrite_bak],1);
if q,
    if system(['cp ',restart_file,' ',rewrite_bak]),
        warning('copy failed')
    end
end

%%% change fuel moisture in restart file %%%
new_moist = w.fmc_g;
new_moist = new_moist.*w_mask;
ncreplace(restart_file,'FMC_G',new_moist)



end % function
