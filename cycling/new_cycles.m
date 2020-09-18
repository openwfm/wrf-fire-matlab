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
fprintf('Adjusting fuels by %f \n',fmc_adjust);
%adjiust fuel gloabally for starters
msk = ones(size(w.xlong));
fmc_change(fmc_adjust,msk);
new_w = insert_analysis(w,ps,tn);
%maybe do this in the fmc_change function?
if cycle == 1
    %back up
    %make more general for more than a single domain simulation
    wi = 'wrfinput_d0';
    domain = input_num('Which wrfinput file to write FMC into? ',1)
    wi = [wi,num2str(domain)];
    wi_bak = [wi,'.bak'];
    cpy_str = sprintf('cp %s %s',wi,wi_bak);
    system(cpy_str)
    ncreplace('wrfinput_d01','TIGN_G',new_w.analysis);
else
    d = dir('wrfrst*');
    for i = 1:length(d)
        fprintf('%d :  %s\n',i,d(i).name)
    end
    r = input_num('Which restart file to use?',1)
    rst = d(r).name;
    rst_bak = [rst,'.bak'];
    cpy_str = sprintf('cp %s %s',rst,rst_bak);
    system(cpy_str);
    %rst='wrfrst_d01_2013-08-13_00:00:00.bak';
    ncreplace(rst,'TIGN_G',new_w.analysis);
end
fprintf('All done. Copy files to directories and restart WRF-SFIRE\n');
fprintf('Do not forget the namelist file needs to be re-linked \n');

end % function
