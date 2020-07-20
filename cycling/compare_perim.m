function compare_perim(w,tign_new)

%w = read_wrfout_tign(f), f - wrfout file
%tign_new = squish(w) - analysis or inital estimate to compare against perimeter


%subset the fire domain to make it easier to work with
red = subset_domain(w);


%% find which fire and data locations to use
sim = input_num('Patch [1] or Camp [2] or Cougar [3] ? ',1,1);
if sim == 1
    det_prefix = '../TIFs/';
    perim = '../PERIMs/patch_perims/';
elseif sim == 2
    det_prefix = '../campTIFs/';
    perim = '../PERIMs/camp_perims/';
else
    det_prefix = '../cougarTIFs/';
    perim = '../PERIMs/cougar_perims/';
end

%% get perim data, arrange and store
fprintf('Reading directory of shape files \n')
a = shape2struct(perim);

%count perims, gather information
%find perimeters from perim file
p_count = 0;

%%%%%%% start 13 for special point %%%%%%%%
%n gives  size of grid to use
n = 100;
for i = 1:length(a)
    %i, a(i)
    if strcmp(a(i).Geometry,'Polygon')
        p_count = p_count + 1;
        %i,a(i)
        
        % get perimeter time
        a(i).p_string = a(i).Name(end-12:end);
        formatIn = 'yyyymmdd HHMM';
        %perim times are local, need to convert to UTC
        %patch fire
        zone_shift = 6;
        % cougar creek , camp fire
        if strcmp(a(i).Name(1:2),'ca') | strcmp(a(i).Name(1:2),'wa')
            %a(i).p_string = a(i).Name(end-12:end);
            zone_shift = 8;
            %formatIn = 'yyyymmdd HHMM';
            %a(i).p_string = a(i).p_string(end-12:end)
        end
        %datenum format as used by TIGN
        a(i).p_time = datenum(a(i).p_string,formatIn)+zone_shift/24;
        
        %set decimate to an  postive integer to use just a subset of points
        %  in perimeter
        fprintf('Perimeter %s: %d points in the perimeter \n',a(i).p_string,length(a(i).Lon));
        if length(a(i).Lat) > 20000
            decimate = 8;
            lats = a(i).Lat(1:decimate:end);
            lons = a(i).Lon(1:decimate:end);
        else
            lats = a(i).Lat;
            lons = a(i).Lon;
        end
        
        %create regularly spaced data
%         dx = (a(i).BoundingBox(2,1)-a(i).BoundingBox(1,1))/n;
%         dy = (a(i).BoundingBox(1,2)-a(i).BoundingBox(2,2))/n;
%         
%         xa = linspace(a(i).BoundingBox(1,1),a(i).BoundingBox(2,1),n);
%         ya = linspace(a(i).BoundingBox(2,2),a(i).BoundingBox(1,2),n);
        
        %find data inside of perimeter
%         [x,y] = meshgrid(xa,ya);
%         x = x(:);
%         y = y(:);
%         [in,on] = inpolygon(x,y,lons,lats);
%         fires = logical(in+on);
%         data = reshape(fires,n,n);
        %make all high confidence fires
%         data = uint8(9.0*data);
%         a(i).data = data;
%         geotransform = [ a(i).BoundingBox(1,1) dx 0  a(i).BoundingBox(2,2) 0 dy];
%         a(i).geotransform = geotransform;
%         %save the file for use in data assimilation
%         %save(a(i).TIF_name,'data','geotransform');
%         %plot results
        plot_on = 1;
        if plot_on
            figure
            contour(red.fxlong,red.fxlat,tign_new,[a(i).p_time a(i).p_time],'b');
            fire_mask = tign_new <= a(i).p_time;
            x_0 = min(red.fxlong(fire_mask));x_1=max(red.fxlong(fire_mask));
            y_0 = min(red.fxlat(fire_mask));y_1=max(red.fxlat(fire_mask));
            hold on
            plot(lons,lats,'r');
            xl_0 = min(lons);xl_1 = max(lons);
            yl_0 = min(lats);yl_1 = max(lats);
            title(a(i).Name);xlabel('Lon'),ylabel('Lat')
            xlim([min(xl_0,x_0) max(xl_1,x_1)]),ylim([min(y_0,yl_0) max(yl_1,y_1)])
            if a(i).p_time > red.end_datenum
                legend('Interpolation forecast','Infrared Perimeter')
            else
            legend('Interpolation','Infrared Perimeter')
            end
            hold off
            %figure,mesh(data)
            
        end %if plot_on
        
        %store perimter structure
        p_struct(p_count).time = a(i).p_time;
        p_struct(p_count).lats = lats;
        p_struct(p_count).lons = lons;
        p_struct(p_count).file = replace(a(i).Name,' ','_');
        p_struct(p_count).Name = a(i).Name;
        
    end
end %for

%save perim_struct.mat a
fprintf('There were %i perimeters found in the data set\n',p_count)

%sort the struct by time first --> last
T = struct2table(p_struct);
sort_table = sortrows(T, 'time');
p_struct = table2struct(sort_table);



end