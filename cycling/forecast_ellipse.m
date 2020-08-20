function e = forecast_ellipse(r,t)

%function e = detection_ellipse(r,t)
% function returns eigenvectors for ellipse containg fire detections
% inputs;
%   r        struct, fire structure from subset_domain
%   t        scalar, t = g(n).time, time of final sateliite data
%output:
%   e        struct with eigenvectors for ellipse contining forecast
%               perimeter

f_mask = r.tign < t;
%fprintf('Forecast are ~= %d \n',sum(f_mask(:)));
f_lons = r.fxlong(f_mask);
f_lats = r.fxlat(f_mask);
forecast = [f_lons,f_lats];

ci = 1.96;
e = ellipse_fit(forecast,ci);
e.area = sum(f_mask(:));
e.ellipse_area = det(e.d);


end


