function fuel_adjust(restart_file)
%function adjust fuel moisture in a  wrf file used for input wrfinput or
%wrfrst

fire_choice = input_num('which fire? Patch: [0], Camp: [1], Cougar: [3]',1)
cycle = input_num('Which cycle? ',0)
load_str = sprintf('p_%d.mat',cycle)
load(load_str);

%nned to beef up what's in red, below
w  = read_wrfout_tign(restart_file);
w.tign_g = p.analysis;
red = subset_domain(w);

%use the wrfout file to find subset of detections

time_bounds(2) = red.max_tign;
end_date = datestr(time_bounds(2));
%time_bounds(2) = 7.354637292824074e+05
time_bounds(1) = red.min_tign;
%time_bounds(1) = time_bounds(2)-3;
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
        fire_mask = g(i).data >=min_con;
        fm = g(i).fxdata >=min_con;
        lons = g(i).xlon(fire_mask);
        lats = g(i).xlat(fire_mask);
    end
    fire_mask = g(i).data >= min_con;
    lon_update = g(i).xlon(fire_mask);
    lat_update = g(i).xlat(fire_mask);
    fm = fm + (g(i).fxdata >=min_con);
    %fprintf('%d detections \n',sum(fm(:)));
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
fm = fm > 0;

forecast_area = red.tign <= g(i).time;
forecast_lons = red.fxlong(forecast_area);
forecast_lats = red.fxlat(forecast_area);

%find perimter of detections
% k = boundary(lons,lats);
% sat_lons = lons(k);
% sat_lats = lats(k);

% [in,on] = inpolygon(forecast_lons,forecast_lats,sat_lons,sat_lats);
% inside = logical(in+on);

mult=1.4;
too_far = mult*(forecast_area-fm);
new2far = ones(size(too_far));
too_far_mask = too_far > 0;
new2far(too_far_mask) = too_far(too_far_mask);

fast_fuel = red.nfuel_cat(too_far_mask);


w_mask = ones(size(w.tign_g));
w_mask(red.ispan,red.jspan) = new2far;

%need to code some way to figure which restart file to use
% if cycle == 1
%     restart = 'wrfinput_d01'
% else
%     restart = 'wrfrst_d01_2013-08-13_00:00:00'
% end

fprintf('%s %s\n','Will write modified time into     ',restart_file)
rewrite_bak=[restart_file,'.bak'];
q=input_num(['1 to copy ',restart_file,' to ',rewrite_bak],1);
if q,
    if system(['cp ',restart_file,' ',rewrite_bak]),
        warning('copy failed')
    end
end

%%% change fuel moisture %%%
new_moist = w.fmc_g;
new_moist = new_moist.*w_mask;
ncreplace(restart_file,'FMC_G',new_moist)



end % function
