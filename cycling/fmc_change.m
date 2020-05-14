function fmc_change(m)

f = 'wrfinput_d01';
%w = read_wrfout_tign(f);
load sm_mask.mat
s = nc2struct('wrfinput_d01',{'FMC_GC'},{})
moist = m*s.fmc_gc;
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
