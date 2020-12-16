function [t,n_diff] = make_tign(ps,alpha,a,b,p_param)
%create palusible fire arrival time from detections
%inputs - 
%    ps = cluster_paths(w,1,grid_size)
%    blending between initial estimate and estimate using paths
%    a,b  - parameters of smooth_up function
%    p_param - parameter of smoothing spline
ps1 = interp_paths(ps,p_param);
%ps1 = ps;
% if ~exist('an_base.mat','file')
%     an = estimate_tign(ps1);
%     save an_base.mat an
% else
%     load an_base.mat
% end
an = estimate_tign(ps1);
%put new estimate in ps1.red
ps1.red.tign = an;
t = squish4(ps1,a,b);

%alpha = 0.4;
t = alpha*t+(1-alpha)*an;
t = smooth_up(ps1.red.fxlat,ps1.red.fxlong,t,a,b);
%t = imgaussfilt(an,1);
figure,mesh(ps.red.fxlong,ps.red.fxlat,t);
%title('Squish4')
hold on
scatter3(ps.points(:,2),ps.points(:,1),ps.points(:,3),'*r');
title('Tign from polygons and paths')
xlabel('Lon'),ylabel('Lat'),zlabel('Datenum')

%compute norms of difference
n_diff = norm(ps.red.tign-t)/norm(ps.red.tign-ps.red.start_datenum);
fprintf('relative error: %f percent\n',100*n_diff)

end %function