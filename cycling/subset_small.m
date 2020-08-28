function r = subset_small(red,n,m)
%interpolates onto smaller grid
r = red;
lat = linspace(r.min_lat,r.max_lat,m);
lon = linspace(r.min_lon,r.max_lon,n);
[r.fxlat,r.fxlong]=meshgrid(lat,lon);

%interpolate dta to smaller grid
r.tign_g = griddata(red.fxlat(:),red.fxlong(:),red.tign_g(:),r.fxlat,r.fxlong);
r.tign = griddata(red.fxlat(:),red.fxlong(:),red.tign(:),r.fxlat,r.fxlong);
r.nfuel_cat = griddata(red.fxlat(:),red.fxlong(:),red.nfuel_cat(:),r.fxlat,r.fxlong);
r.fmc_g = griddata(red.fxlat(:),red.fxlong(:),red.fmc_g(:),r.fxlat,r.fxlong);
r.fhgt = griddata(red.fxlat(:),red.fxlong(:),red.fhgt(:),r.fxlat,r.fxlong);
r.ros = griddata(red.fxlat(:),red.fxlong(:),red.ros(:),r.fxlat,r.fxlong);
r.red = red;




end

