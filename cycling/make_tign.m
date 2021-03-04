function [t,n_diff] = make_tign(ps,alpha,p_param,grid_fraction)
%create palusible fire arrival time from detections
%inputs - 
%    ps = cluster_paths(w,1,grid_size)
%    alpha - blending between initial estimate and estimate made using paths
%    p_param - parameter of smoothing spline
%    grid_fraction - shrink factor fore grid size
ps1 = interp_paths(ps,p_param);
%ps1 = ps;
%%% move to smaller grid
[m,n] = size(ps1.red.fxlong);
%grid_fraction = 0.25
ys = round(grid_fraction*m);
xs = round(grid_fraction*n);
ps1.red = subset_small(ps1.red,ys,xs);
new_points = fixpoints2grid(ps1.red,[ps1.points(:,1), ps1.points(:,2)]);
ps1.idx = new_points(:,1:2);
ps1.grid_pts = new_points(:,3:4);



% if ~exist('an_base.mat','file')
%     an = estimate_tign(ps1);
%     save an_base.mat an
% else
%     load an_base.mat
% end
an = estimate_tign(ps1);
%put new estimate in ps1.red
ps1.red.tign = an;
t = squish4(ps1,1,0);

%alpha = 0.4;
t = alpha*t+(1-alpha)*an;
%t = smooth_up(t,a,b);
%t = imgaussfilt(an,1);
figure,mesh(ps1.red.fxlong,ps1.red.fxlat,t);
%title('Squish4')
hold on
scatter3(ps.points(:,2),ps.points(:,1),ps.points(:,3),'*r');
title('Tign from polygons and paths')
xlabel('Lon'),ylabel('Lat'),zlabel('Datenum')

%resize estimate
if size(t) ~= size(ps.red.fxlong)
    F = scatteredInterpolant(ps1.red.fxlong(:),ps1.red.fxlat(:),t(:));
    t = F(ps.red.fxlong,ps.red.fxlat);
end
%fix the ripple on the top
r_mask = t >=max(t(:))-0.05;
t(r_mask) = max(t(:));

%compute norms of difference
n_diff = norm(ps.red.tign-t)/norm(ps.red.tign-min(ps.red.tign(:)));
fprintf('relative error: %f percent\n',100*n_diff)

end %function
