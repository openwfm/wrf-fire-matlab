function track_struct = track_fmc(wrf)

t=nc2struct(wrf,{'Times'},{});
ts = char(t.times');
[n,m] = size(ts)
for i = 1:n
    s = nc2struct(wrf,{'FMC_G','ROS'},{},i);
    avg_fmc(i) = mean(s.fmc_g(:));
    avg_ros(i) = mean(s.ros(:));
    step_time(i) = datenum(ts(i,:));
end
figure,plot(step_time,avg_fmc)
track_struct.avg_fmc = avg_fmc;
track_struct.step_time = step_time;
track_struct.avg_ros = avg_ros;
end
