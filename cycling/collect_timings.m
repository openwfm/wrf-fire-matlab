function [c,wrf,sfire,write,cores] = collect_timings(test_path)
%collects timings for all wrf runs in a path
%inputs
%   test_path  - string with path to directory whose subdirectories are
%   wrf wruns
%output
%   c  - struct, data structure with timings etc...
r = dir([test_path,'/cyc_i1*']);
n = length(r);
for i = 1:n
   wrf_path = [r(i).folder,'/',r(i).name];
   c(i).name = r(i).name;
   c(i).folder = r(i).folder;
   [c(i).wrf,c(i).sfire,c(i).write] = wrf_timing(wrf_path);
   wrf(i) = c(i).wrf;
   sfire(i) = c(i).sfire;
   write(i) = c(i).write;
   cores(i) = 36;%str2num(c(i).name(end-2:end));
   
end

figure,plot(cores,wrf);
title('Timings for wrf');
xlabel('Cores');ylabel('Time [s]');
figure,plot(cores,sfire);
title('Timings for sfire');
xlabel('Cores');ylabel('Time [s]');

save('timings.mat','c','wrf','sfire','write','cores');


end  % function