function [t,n_diff] = make_tign(ps,alpha)
%create palusible fire arrival time from detections
%inputs - 
%    ps = cluster_paths(w,1,grid_size)
ps1 = interp_paths(ps,0.1);
% if ~exist('an_base.mat','file')
%     an = estimate_tign(ps1);
%     save an_base.mat an
% else
%     load an_base.mat
% end
an = estimate_tign(ps1);
%put new estimate in ps1.red
% ps1.red.tign = an;
% t = squish4(ps1);
t = smooth_up(ps1.red.fxlat,ps1.red.fxlong,an);
%t = imgaussfilt(an,1);
%alpha = 0.4;
t = alpha*t+(1-alpha)*an;
figure,mesh(ps.red.fxlong,ps.red.fxlat,t);
hold on
scatter3(ps.points(:,2),ps.points(:,1),ps.points(:,3),'*r');
title('Tign from polygons and paths')
xlabel('Lon'),ylabel('Lat'),zlabel('Datenum')

%compute norms of difference
n_diff = norm(ps.red.tign-t)/norm(ps.red.tign-ps.red.start_datenum);
fprintf('relative error: %f percent\n',100*n_diff)

end %function