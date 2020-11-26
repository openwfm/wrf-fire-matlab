function r = subset_small(red,n,m)
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
%r.tign_g = griddata(red.fxlat(:),red.fxlong(:),red.tign_g(:),r.fxlat,r.fxlong);
F = scatteredInterpolant(red.fxlat(:),red.fxlong(:),red.tign(:));
r.tign = F(r.fxlat,r.fxlong);
%r.tign = griddata(red.fxlat(:),red.fxlong(:),red.tign(:),r.fxlat,r.fxlong);
F = scatteredInterpolant(red.fxlat(:),red.fxlong(:),red.nfuel_cat(:),'nearest');
r.nfuel_cat = F(r.fxlat,r.fxlong);
%r.nfuel_cat = griddata(red.fxlat(:),red.fxlong(:),red.nfuel_cat(:),r.fxlat,r.fxlong);
F = scatteredInterpolant(red.fxlat(:),red.fxlong(:),red.fmc_g (:),'nearest');
r.fmc_g  = F(r.fxlat,r.fxlong);
%r.fmc_g = griddata(red.fxlat(:),red.fxlong(:),red.fmc_g(:),r.fxlat,r.fxlong);
F = scatteredInterpolant(red.fxlat(:),red.fxlong(:),red.fhgt(:));
r.fhgt = F(r.fxlat,r.fxlong);
%r.fhgt = griddata(red.fxlat(:),red.fxlong(:),red.fhgt(:),r.fxlat,r.fxlong);
F = scatteredInterpolant(red.fxlat(:),red.fxlong(:),red.ros(:));
r.ros = F(r.fxlat,r.fxlong);
%r.ros = griddata(red.fxlat(:),red.fxlong(:),red.ros(:),r.fxlat,r.fxlong);

%trim off NaNs on boundary --> not needed if using scattered interpolant
% r.tign_g = r.tign_g(2:end-1,2:end-1);
% r.tign = r.tign(2:end-1,2:end-1);
% r.nfuel_cat = r.nfuel_cat(2:end-1,2:end-1);
% r.fmc_g = r.fmc_g(2:end-1,2:end-1);
% r.fhgt = r.fhgt(2:end-1,2:end-1);
% r.ros = r.ros(2:end-1,2:end-1);
% r.fxlat = r.fxlat(2:end-1,2:end-1);
% r.fxlong = r.fxlong(2:end-1,2:end-1);
r.red = red;




end

