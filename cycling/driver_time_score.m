%time_score driver
% sets paths to kml files and a wrfout, then runs the time_score function

kml_path = '/bigdisk/james.haley/wrfcycling/wrf-fire/wrfv2_fire/test/TIFs/doc.kml';
wrfout_path = '/bigdisk/james.haley/wrfcycling/wrf-fire/wrfv2_fire/test/cycling_best_ig_single/wrfout_d01_2013-08-13_00:00:00';

time_score(kml_path,wrfout_path)

