function fmc_change(m,msk,f)
%perchentage change to add/subtract
% msk - locations where to add subtract
% f - string, path to a wrfinput or wrfrst file
%f = 'wrfinput_d01';
msk = double(msk);
%msk(msk>0) = 1;

%blur mask a little bit
%msk= imgaussfilt(msk,1/2);
%w = read_wrfout_tign(f);
%load sm_mask.mat
s = nc2struct(f,{'FMC_GC'},{})
fprintf('Fuel levels 1--> 1hr, 2-->10hr 3-->100hr 4-->1000hr 5-->live \n')
fprintf('Standard for now is [3,5]: 100hr, live fuels.\n')
f_time = input_num('Which fuel levels? All = -1',[3,5]);
if f_time < 0
    moist = m*msk + s.fmc_gc;
else
    moist = s.fmc_gc;
    for i = 1:length(f_time)
        %mask area
        moist(:,:,f_time(i)) = m*msk + s.fmc_gc(:,:,f_time(i));
        %small adjustment globally
        moist(:,:,f_time(i)) = moist(:,:,f_time(i))+m/4;
    end
end
%moist2 = moist;
rewrite_bak=[f,'.bak_before_fmc'];
if system(['cp ',f,' ',rewrite_bak])
    fprintf('Error in copy \n')
else 
    fprintf('Copy ok. rewriting FMC_G \n')
    ncreplace(f,'FMC_GC',moist)
    %ncreplace(f,'FMC_G',moist)
end
end