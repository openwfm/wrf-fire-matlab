function new_cycles(f)
% runs new cycling routine as much a spossible
%input:
%    f - string with path to wrfout file

cycle = input_num('Which cycle? ',1);

w = read_wrfout_tign(f);
ps = cluster_paths(w,1);
tn = squish2(ps);
%avg ROS in the forecast and data estimate
[r1,r2] = cone_compare(ps,tn);
%load fuel information from the wrfout
fuels;
figure,plot(fuel(2).fmc_g,fuel(2).ros_fmc_g);
fitter = fit(fuel(2).ros_fmc_g',fuel(2).fmc_g','cubicspline')
%how much to adjust the fmc by
%need to adjust for slope of terrain, right now just use 1/2 of the
%difference
fmc_adjust = 1/2*(fitter(r2)-fitter(r1));
%adjiust fuel gloabally for starters
msk = ones(size(w.xlong));
fmc_change(fmc_adjust,msk);
new_w = insert_analysis(w,ps,tn);
%maybe do this in the fmc_change function?
if cycle == 1
    ncreplace('wrfinput_d01','TIGN_G',new_w.analysis);
else
    fprintf('need to put in selection method for restart file here \n');
    rst='wrfrst_d01_2013-08-13_00:00:00.bak';
    ncreplace(rst,'TIGN_G',new_w.analysis);
end
fprintf('All done. Copy files to directories and restart WRF-SFIRE\n');
fprintf('Do not forget the namelist file needs to be re-linked \n');

end % function
