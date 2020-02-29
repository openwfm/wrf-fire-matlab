function [ score ] = time_score(perim_path,wrfout_path)
%function [ score ] = time_score(perim_kml,wrfout)
% Function computes a score for a fire simulation based on fire arrival
% times of the model compared with satellite observations
% needs the kml2struct function to work
% Inputs:
%       perim - string, path to kml file with perim(s)
%       wrfout - string, path to a wrfout file for a simulation
% Outputs:
%       score - average of the differences in fire arrival times for each
%       pixel in the perimeter

%%%%% read the wrfout file and create interpolant
close all
w = read_wrfout_tign(wrfout_path);
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

%%%%%% convert kml file to a struct
temp_struct = kml2struct(perim_path);

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
       time_str = temp_struct(i).Name(end-13:end);
       perim_date(perim_count) = datetime(time_str,'InputFormat','MM-dd-yyyy HHmm');
       
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

%% compute scores and plot the data

perim_scores = zeros(perim_count,1);
figure(perim_count+1)
mesh(red.fxlong,red.fxlat,(red.tign-red.start_datenum))
%legend('Forecast')
title('Forecast and Perimeters');
xlabel('Lon'),ylabel('Lat'),zlabel('Simulation Time [days]')
for i = 1:perim_count
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
        diff = (z_interp-z)*24;
        perim_scores(i) = mean(abs(diff));
        title_spec = sprintf('Histogram of errors %s',perim_struct(i).Name);
        figure(i),histogram(diff),title(title_spec)
        xlabel('Forecast difference from IR perimeter [hours]')
        ylabel('Number of perimeter points')
        figure(perim_count+1),hold on
        %scatter3(p_lon,p_lat,z_interp,'*')
        scatter3(p_lon,p_lat,(z-red.start_datenum),'.')
        hold off;
        fprintf('%s : Score %f \n', perim_struct(i).Name, perim_scores(i) );
    else
        fprintf('%s perimter after simulation end \n',perim_struct(i).Name);
    end 
    

end
hold off

%%%%%% find the utc times of the perimeters
%% which matlab function for this? 



%%%%%% interpolate perimeter onto the fire mesh
%% use griddata(p_lon,p_lat,z,w.fxlong,w.fxlat) using nearest neighbor
%% z = ones(size(p_lon));

score = mean(perim_scores(perim_scores > 0));

end % function
