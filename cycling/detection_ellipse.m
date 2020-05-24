function e = detection_ellipse(g)
%function e = detection_ellipse(g)
% function returns eigenvectors for ellipse containg fire detections
% inputs;
%   g        struct, detection structure from subset_detections
%output:
%   e        struct with eigenvectors for ellipse contining detections

d = detection_list(g);
fire_idx = find(d);
fires = zeros(length(fire_idx),2);
fire_lons = g(1).xlon(fire_idx);
fire_lats = g(1).xlat(fire_idx);
fires(:,1) = fire_lons;
fires(:,2) = fire_lats;

%scatter plot the detections 
%figure,scatter(fire_lons,fire_lats)

e = ellipse_fit(fires,1.96)



end
