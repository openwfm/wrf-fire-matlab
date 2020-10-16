function g = subset_l2_detections(prefix,p,red,time_bounds,fig)
%same input as load_subset_detections
%reads only active fire pixels, non-detections, water,,etc. are ignored

itime=find(p.time>=time_bounds(1) & p.time<=time_bounds(2));
d=p.file(itime);      % files within the specified time
t=p.time(itime);
fprintf('Selected %i files in the given time bounds, from %i total.\n',...
    length(d),length(p.time))

%granule counter
gc = 0;
for i = 1:length(d)
    file = d{i};
    %only read fire product files
    if strcmp(file(4:5),'14')
        file_str = [prefix,file];

        %fprintf('reading fire product \n')
        if strcmp(file(1),'M')
            try
                fprintf('Reading MODIS FRP data \n')
                v.det = [0 0 1 1 1];
                v.power = hdfread(file_str,'/FP_power');
                v.lon = hdfread(file_str,'/FP_longitude');
                v.lat  = hdfread(file_str,'/FP_latitude');
                v.conf =  hdfread(file_str,'/FP_confidence');
                v.mask = hdfread(file_str,'/fire mask');
            catch
                warning('Read error somewhere')
            end
        else
            fprintf('Reading VIIRS FRP data \n')
            %%fake for the time being
            try
                v.det = [0 0 1 1 1];
                v.power = h5read(file_str,'/FP_power')';
                v.lon = h5read(file_str,'/FP_longitude')';
                v.lat  = h5read(file_str,'/FP_latitude')';
                v.conf =  h5read(file_str,'/FP_confidence')';
                v.mask = h5read(file_str,'/fire mask');
            catch
                warning('read error somewhere')
            end
        end
        %filter data
        xj=find(v.lon > red.min_lon & v.lon < red.max_lon);
        xi=find(v.lat > red.min_lat & v.lat < red.max_lat);
        idx = intersect(xi,xj);
        ax=[red.min_lon red.max_lon red.min_lat red.max_lat];
        if isempty(xi) | isempty(xj)
            fprintf('outside of the domain\n');
        else
            fprintf('inside domain %d \n',gc)
            gc = gc + 1;
            v.lon = v.lon(idx);
            v.lat = v.lat(idx);
            v.power = v.power(idx);
            v.conf = v.conf(idx);
            v.time = t(i);
            v.file = file;
            v.axis=[red.min_lon,red.max_lon,red.min_lat,red.max_lat];
            v.xlon = [];
            v.xlat = [];
            v.fxdata = [];
            %also fake for time being
            v.data = 9*ones(size(v.conf));
            %% put variables into granule struct
            g(gc) = v;
        end
        
    end %
    
    
    
end % for i = ...
fprintf('Detections loaded. %d granules had data inside domain \n',gc);

end % function