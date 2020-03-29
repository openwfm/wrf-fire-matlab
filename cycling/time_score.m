function [ score,diff_var ] = time_score(perim_path,wrfout_path,wrfout_time_step)
%function [ score ] = time_score(perim_kml,wrfout)
% Function computes a score for a fire simulation based on fire arrival
% times of the model compared with satellite observations
% needs the kml2struct function to work
% Inputs:
%       perim - string, path to kml file with perim(s) or to a directory
%          with kmls files or shape file, shape files should be downloaded
%          from geomac as zip, then put into a common directory
%       wrfout - string, path to a wrfout file for a simulation
%       wrfout_time_step - optional, string with the timestep to be read
% Outputs:
%       score - average of the differences in fire arrival times for each
%       pixel in the perimeter

%%%%% read the wrfout file and create interpolant
if nargin > 2
    ts = wrfout_time_step;
    fprintf('Will read wrfout at %s \n',ts)
    w = read_wrfout_tign(wrfout_path,ts);
else
    w = read_wrfout_tign(wrfout_path);
end

%%% compute score using data likelihood
like_score = input_num('Use likelihood approach? No = [0] Yes = [1]',0)
close all

red = subset_domain(w,1);
%% set cone_mask < to use only the cone, no flat part at top
cone_mask = red.tign_g <= max(red.tign_g(:));
flat_top = red.tign_g == max(red.tign_g(:));
%Fw = scatteredInterpolant(w.fxlong(:),w.fxlat(:),w.tign_g(:));
Fr = scatteredInterpolant(red.fxlong(cone_mask),red.fxlat(cone_mask),red.tign(cone_mask));

%%%% testing the interpolation object
%zq = Fr(red.fxlong,red.fxlat);
%figure,mesh(red.fxlong,red.fxlat,red.tign_g),title('org cone')
%figure,mesh(red.fxlong,red.fxlat,zq),title('interpolated new cone')

%%% reading the perimeters
if perim_path(end) == 'l'
    fprintf('Using a kml file for perimeters \n')
    %%%%%% convert kml file to a struct
    temp_struct = kml2struct(perim_path);

    
elseif perim_path(end) == '/'
    
    perim_dat = ['kml';'shp'];
    p_type = input_num('Type of perimeter file to use? (1) kml (2) shp', 1)
    fprintf('Reading %s files in  directory %s \n',perim_dat(p_type,:),perim_path)
    if p_type == 1
        d = dir([perim_path,'*.kml'])
    else
        p_struct = shape2struct(perim_path)
        temp_struct = p_struct;
    end
    
    
end % if perim_path...
    


%%%%%% find the perimeters in the struct, make new struct, find times of
%%%%%% perim in UTC
perim_idx = [];
%perim_time = [];
perim_count = 0;

for i = 1:length(temp_struct)
   % find the perimters, which have the tag 'Polygon'
   if  strcmp(temp_struct(i).Geometry,'Polygon')
       perim_count = perim_count + 1;
       
       perim_idx = [perim_idx,i];
       
       if perim_path(end) == 'l'
          time_str = temp_struct(i).Name(end-13:end);
          perim_date(perim_count) = datetime(time_str,'InputFormat','MM-dd-yyyy HHmm');
       else
          time_str = temp_struct(i).Name(end-12:end);
          perim_date(perim_count) = datetime(time_str,'InputFormat','yyyyMMdd HHmm');
       end
       
       %create new stuct with only perims
       perim_struct(perim_count).Lon = temp_struct(i).Lon;
       perim_struct(perim_count).Lat = temp_struct(i).Lat;
       perim_struct(perim_count).Name = temp_struct(i).Name;
       perim_struct(perim_count).time = datenum(perim_date(perim_count));
       %perim_time = datenum(perim_date);
       %perim_struct(perim_count) = temp_struct(i);
           %perim_struct(perim_count).time = 0;
           % look for nan/inf in data last entry seems to always be a nan....
               %inf_lat = ~isfinite(temp_struct(i).Lat);
               %sum(inf_lat)
       end %if
    end

    %sort the perimeters

    perim_time = datenum(perim_date);
    fprintf('There were %d perimetrs in the kml file \n',perim_count)
    [sorted_times,sort_idx] = sort(perim_date);

    %% compute scores and plot the data

    perim_scores = zeros(perim_count,1);
    perim_vars = zeros(perim_count,1);
    figure(perim_count+1)
    if size(red.fxlong) < 500
        mesh(red.fxlong,red.fxlat,(red.tign-red.start_datenum))
else
    mesh(red.fxlong(1:10:end,1:10:end),red.fxlat(1:10:end,1:10:end),(red.tign(1:10:end,1:10:end)-red.start_datenum))
end
%legend('Forecast')
title('Forecast and Perimeters');
xlabel('Lon'),ylabel('Lat'),zlabel('Simulation Time [days]')
for j = 1:perim_count
    i = sort_idx(j);
    %p_lon = temp_struct(perim_idx(i)).Lon(1:end-1);
    p_lon = perim_struct(i).Lon(1:end-1);
    %p_lat = temp_struct(perim_idx(i)).Lat(1:end-1);
    p_lat = perim_struct(i).Lat(1:end-1);
    
    % only evaluate perims before the end of the simulation
    % scatter plot the perim data
    if perim_struct(i).time <= red.end_datenum - 0.2
        %% ??? perim_time(i)-max(red.max_tign)
        z = perim_struct(i).time*ones(size(p_lon));
        z_interp = Fr(p_lon,p_lat);   
        %%% time difference is perimter time minus forecast time
        %%%   thus early forecast time gives a positive number and 
        %%%   data likelihood would be higher
        diff = (z-z_interp)*24;
        diff = diff(~isnan(diff));
        % test only fire arrival times before simulation end
        cut_top = 0;
        if cut_top
            cone = diff < max(diff)-0.1;
            diff = diff(cone);
        end

        %like_score = 0;
        if like_score > 0
            %make or load spline for data likelihood
            if  ~exist('p_spline','var')
                if exist('p_spline.mat','file')
                    fprintf('Loading spline from mat file \n')
                    load p_spline.mat
                else
                    fprintf('Making and saving spline \n')
                    [p_spline,~,~]= make_spline(72,400);
                    save p_spline.mat p_spline
                end
            else
                fprintf('Using spline in workspace \n')
            end %make spline      
        end %like_score adjustments
        
        
        %%%% output results %%%%
        if like_score == 0
            rel_error = abs(diff(~isnan(diff)))./z(~isnan(diff));
            perim_scores(i) = mean(abs(diff));
            %perim_scores(i) = mean(diff);
            perim_vars(i) = var(diff);
            title_spec = sprintf('Histogram of errors %s',perim_struct(i).Name);
            figure(i),histogram(diff),title(title_spec)
            xlabel('Forecast difference from IR perimeter [hours]')
            ylabel('Number of perimeter points')
            figure(perim_count+1),hold on
            %scatter3(p_lon,p_lat,z_interp,'*')
            scatter3(p_lon,p_lat,(z-red.start_datenum),'.')
            hold off;
            fprintf('%s : Score %f \n', perim_struct(i).Name, perim_scores(i) );
            fprintf('   mean %f var %f \n',mean(diff),var(diff));
            
        else  %%% using likelihood approach
            %using old like
            %stretch=[0.5,10,5,10];
            %dw = ones(size(diff));2
            %likes = like2(dw,diff,stretch);
            likes = p_spline(diff);
            %%exp(likes) give "probabilty"
            %likes = exp(likes);
            perim_scores(i) = mean(likes);
            perim_vars(i) = var(likes);
            title_spec = sprintf('Histogram of likelihoods %s',perim_struct(i).Name);
            figure(i),histogram(likes),title(title_spec)
            xlabel('Perimter point likelihood')
            ylabel('Number of perimeter points')
            figure(perim_count+1),hold on
            %scatter3(p_lon,p_lat,z_interp,'*')
            scatter3(p_lon,p_lat,(z-red.start_datenum),'.')
            hold off;
            fprintf('%s : Score %f \n', perim_struct(i).Name, perim_scores(i) );
            fprintf('   mean %f var %f \n',mean(likes),var(likes));
            
        end
          
    else
        fprintf('%s perimter after simulation end \n',perim_struct(i).Name);
    end 
    

end
hold off

%% use all perims foir the final score ???
% perim_struct.Name
% score_mask = zeros(length(perim_struct),1);
% for k = 1:length(perim_struct)
%     query_string = sprintf('Use %s in score?',perim_struct(k).Name);
%     score_mask(k) = input_num(query_string,1);
% end
% score_mask = logical(score_mask);
    
score_mask = logical([0 0 1 0 0]);
score = perim_scores(score_mask);
diff_var = perim_vars(score_mask);
%score = mean(perim_scores(perim_scores ~= 0));



end % function
