function new_cycles(f)
% runs new cycling routine as much a spossible
%input:
%    f - string with path to wrfout file

cycle = input_num('Which cycle? ',1);

ts = choose_time_step(f);
w = read_wrfout_tign(f,ts);
gs = input_num('What grid spacing?',250)
ps1 = cluster_paths(w,1,gs);
savestr = sprintf('ps_%d.mat',cycle);
save(savestr,'ps1')
%new points
ps = interp_paths(ps1,0.1);
%ps = ps1;

tn = squish4(ps,1,1);
%avg ROS in the forecast and data estimate
[r1,r2,adjr0,outer] = cone_compare(ps,tn);
fprintf('Forecast ROS: %f \n Data ROS: %f \n',r1,r2)
area_diff = outer.a1-outer.a2;
ros_diff = r1-r2;
%load fuel information from the wrfout
fuels;
%figure,plot(fuel(2).fmc_g,fuel(2).ros_fmc_g);
fitter = fit(fuel(2).ros_fmc_g',fuel(2).fmc_g','cubicspline');
%how much to adjust the fmc by
%need to adjust for slope of terrain, right now just use 1/2 of the
%difference
fmc_adjust = 1/2*(fitter(r2)-fitter(r1));
%close all
fprintf('Adjusting fuels by %f percent\n',fmc_adjust);
%adjiust fuel gloabally for starters
new_w = insert_analysis(w,ps,tn);
%maybe do this in the fmc_change function?
if cycle == 1
    %back up
    %make more general for more than a single domain simulation
    wi = 'wrfinput_d0';
    domain = input_num('Which wrfinput_d0x file to write FMC into? x = 1 ',1)
    rst = [wi,num2str(domain)];
    rst_bak = [rst,'.bak'];
    cpy_str = sprintf('cp %s %s',rst,rst_bak);

else
    d = dir('wrfrst*');
    for i = 1:length(d)
        fprintf('%d :  %s\n',i,d(i).name)
    end
    r = input_num('Which restart file to use?',1)
    rst = d(r).name;
    rst_bak = [rst,'.bak'];
    cpy_str = sprintf('cp %s %s',rst,rst_bak);
end

system(cpy_str)
ncreplace(rst,'TIGN_G',new_w.analysis);
%only change fmc content if a bigger fire has a bigger ROS
if ros_diff*area_diff>0
   fmc_change(fmc_adjust,rst);
end
fprintf('Linking new namelist.input file\n')
link_namelist(cycle);
fprintf('All done. Copy files to directories and restart WRF-SFIRE\n');
%fprintf('Do not forget the namelist file needs to be re-linked \n');

end % function
