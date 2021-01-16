function [ scores ] = compare_matches(wrf_single, wrf_cycle, wrf_time_step )
%functio gets comparison scores for two simulations
%inputs:
%  wrf_single, wrf_cycle - strings, paths to wrfouts to be compared
%  wrf_time_step         - string, optional timestep to be read

close all
if nargin > 2
    ts = wrf_time_step;
    scores(1) = match_detections(wrf_single,10,ts);
    scores(2) = match_detections(wrf_cycle,10,ts);
else
    scores(1) = match_detections(wrf_single);
    scores(2) = match_detections(wrf_cycle);
end
    

figure(1), xs = xlim; ys = ylim;
figure(2), xc = xlim; yc = ylim;

x = [xs(:);xc(:)];
y = [ys(:);yc(:)];

x_lim = [min(x(:)) max(x(:))];
y_lim = [min(y(:)) max(y(:))];


figure(1), xlim(x_lim), ylim(y_lim)
figure(2), xlim(x_lim), ylim(y_lim)

end

