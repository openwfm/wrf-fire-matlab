function [i,j,new_lat,new_lon]= fixpt(w,pt,varargin)
% compute the nearest point in the fire mesh to the given coordinates
% w - w = read_wrfout_tign(f)  , struct with fxlong,fxlat or red =
%                                                     subset_domain(w)
%  pt   - point to be fixed to fire mesh pt = [lon,lat]
% 
% i,j  - indices in the grids fxlong,fxlat  or  xlong,xlt
% new_lat,new_lon  -- fxlat(i,j) fxlong(i,j)
% 
% jan mandel, june 2013
% modified by james haley june 2020

lon=w.fxlong;
lat=w.fxlat;

if nargin > 2
    lon = w.xlong;
    lat = w.xlat;
end
min_x=min(lon(:));
max_x=max(lon(:));
min_y=min(lat(:));
max_y=max(lat(:));
[m,n]=size(lon);
[m1,n1]=size(lat);
if(m1 ~= m | n1 ~= n),
    error('inconsistent size of FXLONG and FXLAT')
end
%fprintf('loaded mesh size %i by %i coordinates %g to %g by %g to %g\n',m,n,min_x,max_x,min_y,max_y);


%%%%%%%% compare with what is in w
%find center opf domain
lon_ctr=mean(lon(:));
lat_ctr=mean(lat(:));
%fprintf('the center of the domain is at coordinates %g %g\n',lon_ctr,lat_ctr)
unit_fxlat=6370*2*pi/360;   % one degree latitude in m
unit_fxlon=cos(lat_ctr*2*pi/360)*unit_fxlat; % one degree longitude in m

%fprintf('coordinate units are %g %g m\n',unit_fxlat,unit_fxlon)
    
% find point on the fire mesh
x=pt(2);
y=pt(1);
if (x<min_x | x>max_x | y<min_y | y>max_y),
    error('the point is outside of the domain')
end
% find the nearest point (lon(i,j),lat(i,j)) to (x,y)
d = sqrt((unit_fxlon*(lon-x)).^2 + (unit_fxlat*(lat-y)).^2);
[p,q,minvalue]=find(d==min(d(:)));
for ii=1:length(p)
    i=p(ii);
    j=q(ii);
    new_lon = lon(i,j);
    new_lat = lat(i,j);
    %fprintf('nearest mesh point %i %i at coordinates %12.10g %12.10g distance %g\n',...
    %    i,j,lon(i,j),lat(i,j),d(i,j))
end



end

