function new_cycles(f)
% runs new cycling routine as much a spossible, run from within direcory
% with wrfout file and restart file. Requires fuels.m file
%input:
%    f - string with path to wrfout file
%output:
%   written into restart files

cycle = input_num('Which cycle? ',1);

if strcmp('t',f(end)) %read saved .mat file
    w = read_wrfout_tign(f)
else %read wrfoutfile directly
    ts = choose_time_step(f);
    w = read_wrfout_tign(f,ts);
end

%for handling of branch develop-48 tign_g
t_big = 3.402e+38;
if max(w.tign_g(:)) > t_big
    t_cap = w.itimestep*w.dt;
    msk = w.tign_g >= t_cap;
    w.tign_g(msk) = t_cap;
end

gs = input_num('What grid spacing?',250)
%make the path structure
ps1 = cluster_paths(w,1,gs);
savestr = sprintf('ps_%d.mat',cycle);
save(savestr,'ps1')
%new points interpolated along the paths
ps = interp_paths(ps1,0.9);
%ps = ps1;

%make analysis fire arrival time
tn = squish4(ps,1,1);
%avg ROS in the forecast and data estimate, adjr0 is experimental
%adjustment factor for adjr0 in namelist.fire
[r1,r2,adjr0,outer] = cone_compare(ps,tn);
fprintf('Forecast ROS: %f \n Data ROS: %f \n',r1,r2)
%compare forecast fire area with analysis fire area
area_diff = outer.a1-outer.a2;
ros_diff = r1-r2;
%load fuel information from the wrfout
fuels;
%figure,plot(fuel(2).fmc_g,fuel(2).ros_fmc_g);
%find most common fuel type in burned area
w.nfuel_cat(w.nfuel_cat==14) = NaN;
msk = w.tign_g<max(w.tign_g(:));
common_fuel = mode(w.nfuel_cat(msk));
%find fmc as a function of ROS
fitter = fit(fuel(common_fuel).ros_fmc_g',fuel(common_fuel).fmc_g','cubicspline');
%how much to adjust the fmc by
%need to adjust for slope of terrain, right now just use 1/2 of the
%difference
fmc_adjust = 1/2*(fitter(r2)-fitter(r1));
%close all
fprintf('Adjusting fuels by %f percent\n',fmc_adjust);
%adjiust fuel gloabally for starters

% put analysis fire arrival time struct
new_w = insert_analysis(w,ps,tn);

if cycle == 1
    %restart with wrfinput file
    wi = 'wrfinput_d0';
    domain = input_num('Which wrfinput_d0x file to write FMC into? x = 1 ',1,1)
    rst = [wi,num2str(domain)];
    rst_bak = [rst,'.bak'];
    cpy_str = sprintf('cp %s %s',rst,rst_bak);

else
    %restart with wrfrst file
    d = dir('wrfrst*');
    for i = 1:length(d)
        fprintf('%d :  %s\n',i,d(i).name)
    end
    r = input_num('Which restart file to use?',1)
    rst = d(r).name;
    rst_bak = [rst,'.bak'];
    cpy_str = sprintf('cp %s %s',rst,rst_bak);
end
%copy the file used for the restart
system(cpy_str)
%restart from the analysis
ncreplace(rst,'TIGN_G',new_w.analysis);
%only change fmc content if a bigger fire has a bigger ROS
if ros_diff*area_diff>0
  % fmc_change(fmc_adjust,rst);
end
fprintf('Linking new namelist.input file\n')
link_namelist(cycle);
fprintf('All done. Copy files to directories and restart WRF-SFIRE\n');

end % function
