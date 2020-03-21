function ts = sim_time_scores()

patch_kml = '/bigdisk/james.haley/wrfcycling/wrf-fire/wrfv2_fire/test/TIFs/doc.kml';
time_step = '2013-08-19_00:00:00';
c0 = 'cycle_0/wrfout_d01_2013-08-19_00:00:00';
c1 = 'cycle_1/wrfout_d01_2013-08-19_00:00:00';
c2 = 'cycle_2/wrfout_d01_2013-08-18_00:30:00';
c3 = 'cycle_3/wrfout_d01_2013-08-18_00:30:00';
c4 = 'cycle_4/wrfout_d01_2013-08-18_00:30:00';
c5 = 'wrfout_d01_2013-08-18_00:30:00';

ts = zeros(5,2);
[ts(1,1),ts(1,2)] = time_score(patch_kml,c0);
[ts(2,:),ts(2,2)] = time_score(patch_kml,c1);
[ts(3,:),ts(3,2)] = time_score(patch_kml,c2);
[ts(4,:),ts(4,2)] = time_score(patch_kml,c3);
[ts(5,:),ts(5,2)] = time_score(patch_kml,c4);
[ts(6,:),ts(6,2)] = time_score(patch_kml,c5);
    

end