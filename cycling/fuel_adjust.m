function fuel_adjust(restart_file)
%function adjust fuel moisture in a  wrf file used for input wrfinput or
%wrfrst


w  = read_wrfout_tign(restart_file);
r = subset_domain(w);

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



end % function
