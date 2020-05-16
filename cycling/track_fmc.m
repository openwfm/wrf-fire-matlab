function tracks = track_fmc(wrf)

t=nc2struct(wrf,{'Times'},{});
ts = char(t.times');
for i = 1:length(ts)
    s = nc2struct(wrf,{'FMC_G'},{},i);
    avg_fmc(i) = mean(s.fmc_g(:));
    step_time(i) = datenum(ts(i,:));
end
figure,plot(step_time,avg_fmc)
end
