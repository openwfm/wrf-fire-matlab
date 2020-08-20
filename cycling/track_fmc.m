function track_struct = track_fmc(f)
%compute some trens about fuels, winds, ROS, etc for a wrfout
%inputs
%   f   -   string, path to a wrfout file
%outputs
%   track_struct  - matlab struct with information about the wrfout

t=nc2struct(f,{'Times'},{});
ts = char(t.times');
[n,m] = size(ts)
fuel_cat = 2;
for i = 1:n
    %fprintf('Time : %s \n',ts)
    s = nc2struct(f,{'FMC_G','ROS','NFUEL_CAT','UF','VF','TIGN_G'},{},i);
    fuel_mask = s.nfuel_cat ==2;
    area(i) = sum(sum(s.tign_g<max(s.tign_g(:))));
    avg_fmc(i) = mean(s.fmc_g(fuel_mask));
    avg_ros(i) = mean(s.ros(s.ros>0));
    step_time(i) = datenum(ts(i,:));
    min_fmc(i) = min(s.fmc_g(fuel_mask));
    %compute wind magnitude
    wind_speed = sqrt(s.uf.^2+s.vf.^2);
    %wind(i) = max(max(s.uf(:)),max(s.vf(:)));
    wind(i) = max(wind_speed(:));
    avg_wind(i) = mean(wind_speed(:));
end
%plot reults in days since start
%figure,plot(step_time-step_time(1),avg_fmc)
track_struct.avg_fmc = avg_fmc;
track_struct.step_time = step_time;
track_struct.avg_ros = avg_ros;
track_struct.min_fmc = min_fmc;
track_struct.wind = wind;
track_struct.avg_wind = avg_wind;
track_struct.area = area;
end
