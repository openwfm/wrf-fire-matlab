function e = detection_ellipse(g,r)
%function e = detection_ellipse(g)
% function returns eigenvectors for ellipse containg fire detections
% inputs;
%   g        struct, detection structure from subset_detections
%   r       struct, fire data from subset_domain
%output:
%   e        struct with eigenvectors for ellipse contining detections

d = detection_list(g,r);
fire_idx = find(d);
fires = zeros(length(fire_idx),2);
fire_lons = r.fxlong(fire_idx);
fire_lats = r.fxlat(fire_idx);
fires(:,1) = fire_lons;
fires(:,2) = fire_lats;

%[sat_in,sat_on] = inpolygon(r.fxlong(:),r.fxlat(:),fire_lons,fire_lats);
fire_area = length(fire_idx);
%fprintf('satellite area ~= %d \n',fire_area);

%scatter plot the detections 
% figure,scatter(fire_lons,fire_lats)
% hold on
% shrink = 0.75;
% k = boundary(fire_lons,fire_lats,shrink);
% plot(fire_lons(k),fire_lats(k));

e = ellipse_fit(fires,1.96);
e.area = fire_area;
e.ellipse_area = det(e.d);



end
