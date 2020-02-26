function [ score ] = time_score(perim,wrfout)
%function [ score ] = time_score(perim_kml,wrfout)
% Function computes a score for a fire simulation based on fire arrival
% times of the model compared with satellite observations
% needs the kml2struct function to work
% Inputs:
%       perim - string, path to kml file with perim(s)
%       wrfout - string, path to a wrfout file for a simulation
% Outputs:
%       score - average of the differences in fire arrival times for each
%       pixel in the perimeter

%%%%%% convert kml file to a struct


%%%%%% find the perimeters in the struct


%%%%%% find the utc times of the perimeters



%%%%%% interpolate perimeter onto the fire mesh



%%%%%% compute score for each perimeter









end

