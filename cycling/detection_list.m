function d = detection_list(g)
%d = detection_list(g)
%inputs :
%   g       struct, stelite data from subset detections
%output ;
%   d       nx2  array with lon,lat of detections


min_con = 7;
mask = zeros(size(g(1).data));
for i = 1:length(g)
    mask = mask + g(i).data >= min_con;
end

% plot the fires
%figure,mesh(g(1).xlon,g(1).xlat,mask)

d = mask;
end
    