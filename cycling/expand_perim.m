function ep = expand_perim(red,perim,n)
%function expands the perimeter around a fire
%red = subest_domain(w)
%perim - perim struct from perim2gran function
%n how much blurring, standard deviation in the imgaussfilt function
%find lon,lat within the perimeter
lon = linspace(red.min_lon,red.max_lon,length(perim.lat));
lat = linspace(red.min_lat,red.max_lat,length(perim.lat));
[fxlong,fxlat]=meshgrid(lon,lat);
in = inpolygon(fxlong,fxlat,perim.lon,perim.lat);
%figure,scatter(fxlong(in),fxlat(in))

%blur the mask to expand it
blur = imgaussfilt(double(in),n);
blur = logical(blur);
blur2 = imgaussfilt(double(in),n-1/2);
blur2 = logical(blur2);
%figure,scatter(red.fxlong(new_in),red.fxlat(new_in))
[new_in,new_on] = inpolygon(fxlong,fxlat,fxlong(blur),fxlat(blur));

diff = logical(blur-blur2);
figure,scatter(fxlong(diff),fxlat(diff));
%figure,scatter(red.fxlong(new_on),red.fxlat(new_on));
hold on,scatter(perim.lon,perim.lat)
title('Expanded Perim and Original')

ep = perim;
ep.lon = fxlong(diff);
ep.lat = fxlat(diff);



end %function
