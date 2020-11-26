function outs = l2_detect(red)


[fire_name,save_name,prefix] = fire_choice()
%Longitude = hdfread('/home/jhaley/JPSSdata/MOD03.A2013222.0545.006.2013222112442.hdf', 'MODIS_Swath_Type_GEO', 'Fields', 'Longitude');
%Latitude = hdfread('/home/jhaley/JPSSdata/MOD03.A2013222.0545.006.2013222112442.hdf', 'MODIS_Swath_Type_GEO', 'Fields', 'Latitude');
%fire_mask = hdfread('/home/jhaley/JPSSdata/MOD14.A2013222.0545.006.2015263221706.hdf', '/fire mask', 'Index', {[1  1],[1  1],[2030  1354]});

dhdf=dir([prefix,'*.hdf']);
dh5 = dir([prefix,'*.h5']);
dnc = dir([prefix,'*.nc']);
d=[{dhdf.name},{dh5.name},{dnc.name}];
%d={d.name};

if(isempty(d)), error(['No files found for ',prefix]),end

% order the files in time
nfiles=length(d);
t=zeros(1,nfiles);
for i=1:nfiles
    f{i}=[prefix,d{i}];
    t(i)=rsac2time(d{i});
end

[t,i]=sort(t);
p.file={d{i}};
p.time=t;

fires = [0 0]
gran_count = 1;
for k = 1:nfiles
    if mod(k,2) == 1
        file = p.file{k};
        file2 = p.file{k+1};
        v = readl2data(prefix,file,file2);
        % select fire detection within the domain
        lon_bot = v.lon<red.max_lon;
        lon_top = v.lon>red.min_lon;
        lat_bot = v.lat<red.max_lat;
        lat_top = v.lat>red.min_lat;
        lon_msk = logical(lon_top.*lon_bot);
        lat_msk = logical(lat_top.*lat_bot);
        msk = logical(lat_msk.*lon_msk);
        xj=find(v.lon > red.min_lon & v.lon < red.max_lon);
        xi=find(v.lat > red.min_lat & v.lat < red.max_lat);
        idx = (intersect(xi,xj));
        ax=[red.min_lon red.max_lon red.min_lat red.max_lat];
        if isempty(xi) | isempty(xj)
            fprintf('outside of the domain\n');
        end
        if sum(msk(:)) > 0
        v.data = v.data(idx);
        v.lon = v.lon(msk);
        v.lat = v.lat(msk);
        dets = v.data >= 7;
        g(gran_count) = v;
        v;
        if sum(msk(:)) > 0
            fires = [fires; v.lon(dets) v.lat(dets)];
            %fires = [fires; v.lons(dets) v.lats(dets)];
        end
        gran_count = gran_count + 1;
        end
    end
end

fires = fires(2:end,:);
fprintf('%d fires in %d granules \n',length(fires), length(d))
figure,scatter(fires(:,1),fires(:,2))

outs.fires = fires;
outs.g = g;

end

