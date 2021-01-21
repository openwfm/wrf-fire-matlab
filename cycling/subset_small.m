function r = subset_small(red,n,m,full_set)
%interpolates onto smaller grid
% red = subset_domain(w)
% n,m new size of matrices
r = red;
% d_lat = (red.max_lat-red.min_lat)/m;
% d_lon = (red.max_lon-red.min_lon)/n;
% lat = linspace(r.min_lat-d_lat,r.max_lat+d_lat,m+2);
% lon = linspace(r.min_lon-d_lon,r.max_lon+d_lon,n+2);
lon = linspace(r.min_lon,r.max_lon,n);
lat = linspace(r.min_lat,r.max_lat,m);
[r.fxlat,r.fxlong]=meshgrid(lat,lon);

%interpolate dta to smaller grid
F = scatteredInterpolant(red.fxlat(:),red.fxlong(:),red.tign_g(:));
r.tign_g = F(r.fxlat,r.fxlong);

F = scatteredInterpolant(red.fxlat(:),red.fxlong(:),red.tign(:));
r.tign = F(r.fxlat,r.fxlong);

if ~exist('full_set','var')
    full_set = input_num('Interpolate full set of values? Yes = 1',0,1);
end
if full_set == 1
    F = scatteredInterpolant(red.fxlat(:),red.fxlong(:),red.nfuel_cat(:),'nearest');
    r.nfuel_cat = F(r.fxlat,r.fxlong);
    
    F = scatteredInterpolant(red.fxlat(:),red.fxlong(:),red.fmc_g (:),'nearest');
    r.fmc_g  = F(r.fxlat,r.fxlong);
    
    F = scatteredInterpolant(red.fxlat(:),red.fxlong(:),red.fhgt(:));
    r.fhgt = F(r.fxlat,r.fxlong);
    
    F = scatteredInterpolant(red.fxlat(:),red.fxlong(:),red.ros(:));
    r.ros = F(r.fxlat,r.fxlong);
end

r.red = red;




end

