function [wf,wa,ws] = p2w(p,w)
% function takes in a wrfout file, stored at matlab .mat and returns the
% structs using the forecast, analysis, and spinup from p =
% detect_fit_level2

wf = w;
wa = w;
wa.tign_g = p.analysis;
ws = w;
ws.tign_g = p.spinup;


end