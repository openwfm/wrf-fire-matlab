function score = evaluate_sim( wrfout_s, wrfout_c, wrf_time_step )
%score = evaluate_sim( wrfout, perim )
% inputs:
%   wrfout_s - wrfout from a simulation w/o cycling or spinup or ROS adjustment, string with file path
%   wrfout_c - wrfout from a simulation with cycling or spinup or ROS adjustment, string with file path
%   wrf_time_step - optional string with time step to be read
% output:
%   score - number which gives how good the simulation matches the perimter
%   data
% note: needs the function kml2mat to run

%to do
%
% scatter detections first
% use alpha on detections?
% use wrfout with earliest time comparable to perimeter time by searching
% through the time steps of the wrfouts and comparing with the perimeter
% time. input will then be directory and domain level instead of wrfout
% name

close all

min_fire_confidence = 7;

%exp_num = input_num('Compare No ROS adjust vs ROS adjust = [1] Cycle vs No cycle = [2] Two diff cycles = [3]',1,0);
exp_num = input_num('Compare old method vs new method = [1] Cycle vs No cycle = [2] Two diff cycles = [3]',1,0);
if exp_num == 3
   cycle_s = input_num('Cycle number for wrfout_s?',0);
   cycle_c = input_num('Cycle number for wrfout_c?',1);
end
%read wrfout section
if nargin > 2
    w_c =read_wrfout_tign(wrfout_c,wrf_time_step);
    red_c = subset_domain(w_c,1);
    clear('w_c')
    w_s =read_wrfout_tign(wrfout_s,wrf_time_step);
    red_s = subset_domain(w_s,1);
    clear('w_s');
else
    w_c =read_wrfout_tign(wrfout_c);
    red_c = subset_domain(w_c);
    clear('w_c');
    w_s =read_wrfout_tign(wrfout_s);
    red_s = subset_domain(w_s);
    clear('w_s');
end



%collect information about simulations
max_tign = max(red_s.max_tign,red_c.max_tign);
min_tign = min(red_s.min_tign,red_c.min_tign);

% need to plot detections too, this will be moved elsewhere
sim = input_num('Patch [1] or Camp [2] or Cougar [3] ? ',1,1);
if sim == 1
    det_prefix = '../TIFs/';
    perim = '../PERIMs/patch_perims/';
elseif sim == 2
    det_prefix = '../campTIFs/';
    perim = '../PERIMs/camp_perims/';
else
    det_prefix = '../cougarTIFs/';
    perim = '../PERIMs/cougar_perims/';
end

det_list=sort_rsac_files(det_prefix);
fig.fig_map=0;
fig.fig_3d=0;
fig.fig_interp=0;


if perim(end) == 'l'
    %read kml section, *.kml file
    a = kml2struct(perim);
else
    % or
    fprintf('Reading directory of shape files \n')
    a = shape2struct(perim);
end

%find perimeters from perim file
p_count = 0;

%%%%%%% start 13 for special point %%%%%%%%
%n gives  size of grid to use
n = 100;
for i = 1:length(a)
    %i, a(i)
    if strcmp(a(i).Geometry,'Polygon')
        p_count = p_count + 1;
        %i,a(i)
        
        % get perimeter time
        a(i).p_string = a(i).Name(end-12:end);
        
%        only for kml file
%         if a(i).p_string(1) ~= '1'
%             a(i).p_string(1) = '0';
%         end
%         formatIn = 'mm-dd-yyyy HHMM';
        formatIn = 'yyyymmdd HHMM';
        %perim times are local, need to convert to UTC
        zone_shift = 6;
        if strcmp(a(i).Name(1:2),'ca') | strcmp(a(i).Name(1:2),'wa')
            %a(i).p_string = a(i).Name(end-12:end);
            zone_shift = 8;
            %formatIn = 'yyyymmdd HHMM';
            %a(i).p_string = a(i).p_string(end-12:end)
        end
        a(i).p_time = datenum(a(i).p_string,formatIn)+zone_shift/24;
        
        %set decimate to an  postive integer to use just a subset of points
        %  in perimeter
        fprintf('Perimeter %s: %d points in the perimeter \n',a(i).p_string,length(a(i).Lon));
        if length(a(i).Lat) > 500
            decimate = 10;
            lats = a(i).Lat(1:decimate:end);
            lons = a(i).Lon(1:decimate:end);
        else
            lats = a(i).Lat;
            lons = a(i).Lon;
        end
        
        %create regularly spaced data
        dx = (a(i).BoundingBox(2,1)-a(i).BoundingBox(1,1))/n;
        dy = (a(i).BoundingBox(1,2)-a(i).BoundingBox(2,2))/n;
        
        xa = linspace(a(i).BoundingBox(1,1),a(i).BoundingBox(2,1),n);
        ya = linspace(a(i).BoundingBox(2,2),a(i).BoundingBox(1,2),n);
        
        %find data inside of perimeter
        [x,y] = meshgrid(xa,ya);
        x = x(:);
        y = y(:);
        [in,on] = inpolygon(x,y,lons,lats);
        fires = logical(in+on);
        data = reshape(fires,n,n);
        %make all high confidence fires
        data = uint8(9.0*data);
        a(i).data = data;
        geotransform = [ a(i).BoundingBox(1,1) dx 0  a(i).BoundingBox(2,2) 0 dy];
        a(i).geotransform = geotransform;
        %save the file for use in data assimilation
        %save(a(i).TIF_name,'data','geotransform');
        %plot results
        plot_on = 0;
        if plot_on
            figure
            hold on
            plot(lons,lats);
            scatter(x,y);
            scatter(x(fires),y(fires),'*','r');
            title(a(i).Name);
            %plot(lons,lats)
            hold off
            %figure,mesh(data)
        end %if plot_on
        
        %store perimter structure
        p_struct(p_count).time = a(i).p_time;
        p_struct(p_count).lats = lats;
        p_struct(p_count).lons = lons;
        p_struct(p_count).file = replace(a(i).Name,' ','_');
        p_struct(p_count).Name = a(i).Name;
        
    end
end %for

%save perim_struct.mat a
fprintf('There were %i perimeters found in the data set\n',p_count)

%sort the struct by time first --> last
T = struct2table(p_struct);
sort_table = sortrows(T, 'time');
p_struct = table2struct(sort_table);


%compare perimiters with the wrfout
% need to filter out perimeters which are before first ignition time, score
% would be zero


for i = 1:p_count
    if p_struct(i).time <= max_tign && p_struct(i).time >= min_tign
        fprintf('Comparing wrfout files with the infrared perimeter \n')
        
        sim_fires_s = red_s.tign <= p_struct(i).time;
        sim_fires_c = red_c.tign <= p_struct(i).time;
        
        max_lat_s = max(red_s.fxlat(sim_fires_s));
        min_lat_s = min(red_s.fxlat(sim_fires_s));
        max_lat_c = max(red_c.fxlat(sim_fires_c));
        min_lat_c = min(red_c.fxlat(sim_fires_c));
        max_lon_s = max(red_s.fxlong(sim_fires_s));
        min_lon_s = min(red_s.fxlong(sim_fires_s));
        max_lon_c = max(red_c.fxlong(sim_fires_c));
        min_lon_c = min(red_c.fxlong(sim_fires_c));
        
        %put both tign and perim on the same grid
        max_lat = max([max_lat_s,max_lat_c,max(p_struct(i).lats(:))])+0.01;
        min_lat = min([min_lat_s,min_lat_c,min(p_struct(i).lats(:))])-0.01;
        max_lon = max([max_lon_s,max_lon_c,max(p_struct(i).lons(:))])+0.01;
        min_lon = min([min_lon_s,min_lon_c,min(p_struct(i).lons(:))])-0.01;
        
        mesh_dim = 200;
        x_new = linspace(min_lon,max_lon,mesh_dim);
        y_new = linspace(min_lat,max_lat,mesh_dim);
        %scatter(x_new,y_new)
        [x_grid,y_grid] = meshgrid(x_new,y_new);
        x_vect = x_grid(:);
        y_vect = y_grid(:);
        tign_grid_s = griddata(red_s.fxlong,red_s.fxlat,red_s.tign,x_grid,y_grid);
        tign_grid_c = griddata(red_c.fxlong,red_c.fxlat,red_c.tign,x_grid,y_grid);
        %figure,mesh(x_grid,y_grid,tign_grid)
        
        %mask of fire areas
        sim_fires_s = tign_grid_s <= p_struct(i).time;
        sim_fires_c = tign_grid_c <= p_struct(i).time;
        
        
        %compute areas of simulations
        area_sim_s = sum(sim_fires_s(:));
        x_fires_sim_s = x_grid(sim_fires_s);
        y_fires_sim_s = y_grid(sim_fires_s);
        area_sim_c = sum(sim_fires_c(:));
        x_fires_sim_c = x_grid(sim_fires_c);
        y_fires_sim_c = y_grid(sim_fires_c);
        shrink = 0.95;
        %find boundaries
        sim_boundary_s = boundary(x_fires_sim_s,y_fires_sim_s,shrink);
        sim_x_s = x_fires_sim_s(sim_boundary_s);
        sim_y_s = y_fires_sim_s(sim_boundary_s);
        sim_boundary_c = boundary(x_fires_sim_c,y_fires_sim_c,shrink);
        sim_x_c = x_fires_sim_c(sim_boundary_c);
        sim_y_c = y_fires_sim_c(sim_boundary_c);
        
        %mesh(sim_fires)
        
        % put perim on the grid in use...
        [in_p,on_p] = inpolygon(x_grid(:),y_grid(:),p_struct(i).lons(:),p_struct(i).lats(:));
        perim_fires = logical(in_p+on_p);
        area_perim = sum(perim_fires(:));
        x_fires_perim = x_grid(perim_fires);
        y_fires_perim = y_grid(perim_fires);
        %shrink = 0.75; already declared
        perim_boundary = boundary(x_fires_perim,y_fires_perim,shrink);
        perim_x = x_fires_perim(perim_boundary);
        perim_y = y_fires_perim(perim_boundary);
        
        perim_fires_m = reshape(perim_fires,mesh_dim,mesh_dim);
        
        %calculate score
        diff_s = abs(perim_fires_m - sim_fires_s);
        diff_c = abs(perim_fires_m - sim_fires_c);
        %figure, mesh(diff);
        area_diff_s = sum(diff_s(:));
        area_diff_c = sum(diff_c(:));
        
        score_c(i) = (area_sim_c+area_perim-area_diff_c)/(area_sim_c+area_perim);
        score_s(i) = (area_sim_s+area_perim-area_diff_s)/(area_sim_s+area_perim);
        
        %collect detection data here
        if i == 1
            time_bounds(1) = min_tign;
        else 
            time_bounds(1) = p_struct(i-1).time;
        end
        time_bounds(2) = p_struct(i).time; % will use perimeter time
        g_c = load_subset_detections(det_prefix,det_list,red_c,time_bounds,fig);
        %g_s = load_subset_detections(det_prefix,det_list,red_s,time_bounds,fig);
        % need to loop throught the g_* files and collect all detections

        for j = 1:length(g_c)
            if j == 1
                det_g_c = g_c(j).data >= min_fire_confidence;
                scatter_lon = g_c(j).xlon(det_g_c);
                scatter_lat = g_c(j).xlat(det_g_c);
            end
            det_g_c = g_c(j).data >= min_fire_confidence;
            lon_g_c = g_c(j).xlon(det_g_c);
            lat_g_c = g_c(j).xlat(det_g_c);
            %scatter_det = [scatter_det det_g_c(:)];
            if size(lon_g_c) > 100
                %decimate = uint8(log(size(lon_g_c)));
                decimate = 2;
            end
            % if more than 100 detections, plot just a subset
            if sum(det_g_c(:)) > 0
                scatter_lon = [scatter_lon(:);lon_g_c(1:decimate:end)];
                scatter_lat = [scatter_lat(:);lat_g_c(1:decimate:end)];
            end
        end        
        %plot perims w/o cycling
        figure(2*i);
        hold on
        %scatter detections here
        alpha(0.5);
        scatter(scatter_lon,scatter_lat,'m','*');        alpha(1.0);
        plot(sim_x_s,sim_y_s,'b');
        plot(perim_x ,perim_y,'r');
        %xticks(-112.72:0.02:-112.58)
        if exp_num == 2
            title_str = ('Perimeter Observation and Forecast without cycling');
            legend({'Satellite Fire Detections','Forecast without cycling','Infrared perimeter'});
            save_str = [p_struct(i).file '_s']; %% _s is for single run
        elseif exp_num == 1
              title_str = ('Perimeter Observation and Forecast using old method');
              legend({'Satellite Fire Detections','Forecast from old method','Infrared perimeter'});
              save_str = [p_struct(i).file '_old'];
%             title_str = ('Perimeter Observation and Forecast without spin-up');
%             legend({'Satellite Fire Detections','Forecast without spin-up','Infrared perimeter'});
%             save_str = [p_struct(i).file '_No_spinup'];
%             title_str = ('Perimeter Observation and Forecast without ROS  adjust');
%             legend({'Satellite Fire Detections','Forecast without ROS adjust','Infrared perimeter'});
%             save_str = [p_struct(i).file '_No_ROS_adj'];
        else %spinup_compare == 3
            title_str = sprintf('Perimeter Observation and Forecast for Cycle %d ',cycle_s);
            %title_str = ('Perimeter Observation and Forecast without spin-up');
            legend({'Satellite Fire Detections','Forecast','Infrared perimeter'});
            save_str = sprintf('%s_cycle_%d',p_struct(i).file,cycle_s);
        end
        
        %legend({'Forecast without cycling','Infrared perimeter'});
        xlabel('Lon')
        ylabel('Lat')
        xlim([min_lon max_lon])
        ylim([min_lat max_lat])
        
        score_str = sprintf('Score = %f',score_s(i));
        title({title_str,p_struct(i).Name,score_str});
%         if spinup_compare ~= 1
%             save_str = [p_struct(i).file '_s'];
%         else
%             save_str = [p_struct(i).file '_No_spinup'];
%         end
        
        savefig(save_str);
        saveas(gcf,[save_str '.png']);
        ax(i,1) = gca;
        hold off
        %plot perims w/ cycling
        figure(2*i+1);
        hold on
        %scatter detections here
        alpha(0.5);
        scatter(scatter_lon,scatter_lat,'m','*');
        alpha(1.0);
        plot(sim_x_c,sim_y_c,'b');
        plot(perim_x ,perim_y,'r');
        %xticks(-112.72:0.02:-112.58)
        if exp_num == 2
            legend({'Satellite Fire Detections','Forecast with cycling','Infrared perimeter'});
            title_str = ('Perimeter Observation and Forecast with cycling');
            save_str = [p_struct(i).file '_c'];
        elseif exp_num == 1
              title_str = ('Perimeter Observation and Forecast using new method');
              legend({'Satellite Fire Detections','Forecast from new method','Infrared perimeter'});
              save_str = [p_struct(i).file '_new'];
%             legend({'Satellite Fire Detections','Forecast using spin-up','Infrared perimeter'});
%             title_str = ('Perimeter Observation and Forecast using spin-up');
%             save_str = [p_struct(i).file '_with_spinup'];
%             legend({'Satellite Fire Detections','Forecast using ROS adjust','Infrared perimeter'});
%             title_str = ('Perimeter Observation and Forecast using ROS adjust');
%             save_str = [p_struct(i).file '_with_ROS_adj'];
        else %spinup_compare == 3
            title_str = sprintf('Perimeter Observation and Forecast for Cycle %d ',cycle_c);
            %title_str = ('Perimeter Observation and Forecast without spin-up');
            legend({'Satellite Fire Detections','Forecast','Infrared perimeter'});
            save_str = sprintf('%s_cycle_%d',p_struct(i).file,cycle_c);
        end
        
        %legend({'Forecast with cycling','Infrared perimeter'});
        xlabel('Lon')
        ylabel('Lat')
        xlim([min_lon max_lon])
        ylim([min_lat max_lat])
%         if spinup_compare ~= 1
%             save_str = [p_struct(i).file '_c'];
%         else
%             save_str = [p_struct(i).file '_with_spinup'];
%         end
            
        score_str = sprintf('Score = %f',score_c(i));
        title({title_str,p_struct(i).Name,score_str});
        savefig(save_str);
        saveas(gcf,[save_str '.png']);
        ax(i,2) = gca;
        hold off
        
    else
        if p_struct(i).time < min_tign
            fprintf('Minimum tign after perimeter time, score is zero. Not plotting. \n');
            score_s(i) = 1e-16;
            score_c(i) = 1e-16;
            %             sim_x_s = [];
            %             sim_y_s = [];
            %             sim_x_c = [];
            %             sim_y_c = [];
        else
            fprintf('Perimiter time after simulation end. No score assigned. \n')
        end
    end
    
    
end


%do one plot with everything
% [m,n] =size(ax)
% fnew = figure;
% for j = 1:m
%     ax_copy(j,1) = copyobj(ax(j,1),fnew);
%     ax_copy(j,2) = copyobj(ax(j,2),fnew);
%     subplot(m,n,2*j-1,ax_copy(j,1))
%     subplot(m,n,2*j,ax_copy(j,2))
%     
%     
% end %for j

whos
%m_s = score_s > 0;
m_c = score_c > 0.02;

mean_score_s = mean(score_s(m_c));
mean_score_c = mean(score_c(m_c));


score.score_s = score_s;
score.score_c = score_c;
score.mean_score_s = mean(score_s(m_c));
score.mean_score_c = mean(score_c(m_c));
score.file = p_struct.file;
end

