%% DOCUMENT TITLE
% INTRODUCTORY TEXT
%%
function a = area_ellipse(wrfout)
% function a = area_ellipse(wrfout)
% function for comparing area and orientations
%       of stellite data and model forecasts
% input:
%   wrfout     -  string, path to  wrfout file from wrf
% ouput:
%   a             struct with who knows what in it

[fire_name,save_name,prefix] = fire_choice()

w = read_wrfout_tign(wrfout);
red = subset_domain(w);

% figures
fig.fig_map=0;
fig.fig_3d=0;
fig.fig_interp=0;


time_bounds = [red.start_datenum red.end_datenum];
p = sort_rsac_files(prefix);
if ~exist('g_ellipse.mat','file')
    g = load_subset_detections(prefix,p,red,time_bounds,fig);
    save g_ellipse.mat g
else
    load_new = input_num('re-load detections from scratch? ',1);
    if load_new
        g = load_subset_detections(prefix,p,red,time_bounds,fig);
        save g_ellipse.mat g
    else
        load g_ellipse.mat
    end
end

t = g(end).time;

%fite ellipse, plot results
figure(1)
e_data = detection_ellipse(g,red);
hold on
e_forecast = forecast_ellipse(red,t);
contour(red.fxlong,red.fxlat,red.tign,'k')

a.e_data = e_data;
a.e_forecast = e_forecast;

eig_data = diag(a.e_data.d);
eig_forecast = diag(a.e_forecast.d);
eig_ratio = eig_forecast./eig_data;

fprintf('Forecast Area / Sat Area  = %f \n',e_forecast.area/e_data.area);
fprintf('Forecast det  / Sat det   = %f \n',e_forecast.ellipse_area/e_data.ellipse_area);
fprintf('Forecast eigenvalues / Sat eignevalues %f   %f   \n',eig_ratio(1),eig_ratio(2));


end



