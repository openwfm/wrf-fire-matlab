function track_struct = track_fmc(wrf)

t=nc2struct(wrf,{'Times'},{});
ts = char(t.times');
[n,m] = size(ts)
fuel_cat = 2;
for i = 1:n
    s = nc2struct(wrf,{'FMC_G','ROS','NFUEL_CAT','UF','VF'},{},i);
    fuel_mask = s.nfuel_cat ==2;
    avg_fmc(i) = mean(s.fmc_g(:));
    avg_ros(i) = mean(s.ros(:));
    step_time(i) = datenum(ts(i,:));
    min_fmc(i) = min(s.fmc_g(fuel_mask));
    wind(i) = max(max(s.uf(:)),max(s.vf(:)));
end
figure,plot(step_time,avg_fmc)
track_struct.avg_fmc = avg_fmc;
track_struct.step_time = step_time;
track_struct.avg_ros = avg_ros;
track_struct.min_fmc = min_fmc;
track_struct.wind = wind;
end
