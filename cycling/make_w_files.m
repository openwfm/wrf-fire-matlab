function wrfs = make_w_files(varargin)
%make matlab file from wrfouts

if nargin > 0
    cwd = varargin{1};
else
    cwd = pwd;
end

d = dir(cwd);
% d(1).name
% d(2).name
%all directories have '.' and '..' as dir(1),dir(2)??
for i = 3:length(d)
    if d(i).isdir
        dir_name = ([d(i).name,'/']);
        dir_fold = ([d(i).folder,'/']);
        dir_path = ([dir_fold,dir_name]);
        %cd(dir_path)
        fprintf('calling make_w_files on %s \n',dir_path)
        make_w_files(dir_path)
        %cd('../')
        %fprintf('moved to %s \n',pwd)
        
        %%%strcmp can avoid this...
    elseif length(d(i).name) > 6
        if d(i).name(1:6) == 'wrfout'
            
            %dir_name = ([d(i).name,'/']);
            wrf_fold = ([d(i).folder,'/']);
            wrf_path = ([wrf_fold,d(i).name]);
            
            w = read_wrfout_tign(wrf_path);
            save_str = ([wrf_path,'.mat']);
            fprintf('Working in %s \n', cwd)
            fprintf('Saving %s \n',save_str)
            save(save_str,'w');
        end
        %
    end

end