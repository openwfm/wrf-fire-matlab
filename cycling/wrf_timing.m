function [wrft,sfire,write] = wrf_timing(wrf_path)
% function [wrf,sfire,write] = wrf_timing(wrf_path)
% function computes timing from rsl files in directory wrf_path

% count the rsl.out files
outs = dir([wrf_path,'/rsl*err*']);
l = length(outs);

wrt_total = 0;
wrf_total = 0;
sfr_total = 0;

%loop through the rls.out files, sum the timings

%loop from 1:1 since rsl.error.oooo has all the timings

for i = 1:1
    f = [outs(i).folder,'/',outs(i).name];
    fprintf('Processing file %s \n',f)
    
    wrt_file = [outs(i).folder,'/writ.txt']; %'writ.txt';
    wrf_file = [outs(i).folder,'/wrf.txt'];%'wrf.txt';
    sfr_file = [outs(i).folder,'/wfr.txt'];%'sfr.txt';
    
    %store the different timings
    write_str = sprintf('grep "Timing for Writing" %s > %s',f,wrt_file);
    system(write_str);
    wrf_str = sprintf('grep "Timing for main" %s > %s',f,wrf_file);
    system(wrf_str);
    sfr_str = sprintf('grep "Timing for sfire" %s > %s',f,sfr_file);
    system(sfr_str);
    
    wrt = fopen(wrt_file,'r');
    wrt_spec = '%s %s %s %s %s %s %s %f %s %s';
    C_wrt = textscan(wrt,wrt_spec);
    wrt_time = C_wrt{8};
    fclose(wrt);
    wrt_sum(i) = sum(wrt_time);
    wrt_total = wrt_sum(i) + wrt_total;
    fprintf('Total writing time %s: %f\n',f,wrt_sum(i));
    
    wrf = fopen(wrf_file,'r');
    wrf_spec = '%s %s %s %s %s %s %s %s %f %s %s';
    C_wrf = textscan(wrt,wrf_spec);
    wrf_time = C_wrf{9};
    fclose(wrf);
    wrf_sum(i) = sum(wrf_time);
    wrf_total = wrf_sum(i) + wrf_total;
    fprintf('Total main time %s: %f\n',f,wrf_sum(i));
    
    sfr = fopen(sfr_file,'r');
    sfr_spec = '%s %s %s %f %s %s';
    C_sfr = textscan(sfr,sfr_spec);
    sfr_time = C_sfr{4};
    fclose(sfr);
    sfr_sum(i) = sum(sfr_time);
    sfr_total = sfr_sum(i) + sfr_total;
    fprintf('Total sfire time %s: %f\n',f,sfr_sum(i));


end % for i

wrft = wrf_total;
sfire = sfr_total;
write = wrt_total;

end % function