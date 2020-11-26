function [wf,wa,ws] = p2w(p,w)
% function takes in a wrfout file, stored at matlab .mat and returns the
% structs using the forecast, analysis, and spinup from p =
% detect_fit_level2

wf = w;
save w_forecast.mat w;

w.tign_g = p.analysis;
wa = w;
save w_analysis.mat w;

w.tign_g = p.spinup;
ws = w;
save w_spinup.mat w;



end