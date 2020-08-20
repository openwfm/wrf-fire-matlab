function d = detection_list(g,r)
%d = detection_list(g)
%inputs :
%   g       struct, satellite data from subset_detections
%   r       struct, fire data from subset_domain
%output ;
%   d       nx2  array with lon,lat of detections


min_con = 7;
mask = zeros(size(g(1).fxdata));
for i = 1:length(g)
    mask = mask + double(g(i).fxdata >= min_con);
end

% plot the fires
%figure,mesh(r.fxlong,r.fxlat,mask)

d = mask;
end
    