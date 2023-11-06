
function ts=ts_at(lon,lat,xlon,xlat,v)
% interpolate time series of a variable from grid to a point
% in:
%   lon longitude of the location to interpolate to
%   lat latitude of the location to interpolate to
%   xlon 2D array of longitude coordinates at grid points
%   xlat 2D array of latitude coordinates of grid points
%   v 3D array values. The last index is time.
% out:
%  ts time series of v interpolated to (lon,lat), 1D array

% Method: Interpolate on a small square subgrid around the point (lon,lat). 
% It is NOT assumed that the grid is equidistant, for example xlong xlat
% can be coordinates in one projection while the grid is equidistant in
% another projection. The interpolation is exact if v is linear function
% of xlong and xlat.

function ts_at_test
% generate xlong xlat as uniform with random perturbations
% define v=a*xlong + b*xlat for some a,b
% verify that ts_at returns a*lon + b*lat for several random (lon,lat)
