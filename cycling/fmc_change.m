function fmc_change(m)

f = 'wrfinput_d04';
%w = read_wrfout_tign(f);
%load sm_mask.mat
s = nc2struct(f,{'FMC_GC'},{})
f_time = input_num('Which fuel levels? All = -1',-1);
if f_time < 0
    moist = m*s.fmc_gc;
else
    moist = s.fmc_gc;
    for i = 1:length(f_time)
        moist(:,:,f_time(i)) = m*s.fmc_gc(:,:,f_time(i));
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
