function new_pts = fixpoints2grid(r,pts)
%function moves scattered points to closest locations on fire grid
% inputs:
%    r,    red,w   --   struct comtaining fxlong,fxlat from a wrfout or
%    wrfinput
%    pts  nx2 matrix with rows [lat lon]
% output:
%    new_pts  - indices i,j for closest locations on the fire mesh
%    fxlong,fxlat
new_pts = zeros(length(pts),4);
for i = 1:length(pts)
    [new_pts(i,1),new_pts(i,2),new_pts(i,3),new_pts(i,4)] = fixpt(r,pts(i,1:2));
end
    