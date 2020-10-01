function p_struct = perim2gran(n,perim)
%reads in a directory of perimieters and makes them into "satellite granules"
%n - number of points from each perimeter to use
%perim - path to the perimeter data


%%% reading the perimeters
if perim(end) == 'l'
    fprintf('Using a kml file for perimeters \n')
    %%%%%% convert kml file to a struct
    temp_struct = kml2struct(perim);
    
    
elseif perim(end) == '/'
    
    perim_dat = ['kml';'shp'];
    p_type = input_num('Type of perimeter file to use? (1) kml (2) shp', 2,1);
    fprintf('Reading %s files in  directory %s \n',perim_dat(p_type,:),perim);
    if p_type == 1
        d = dir([perim,'*.kml']);
    else
        a = shape2struct(perim);
    end
end % if perim_path...

%count perims, gather information
%find perimeters from perim file
p_count = 0;

%n gives  number of points in the perimeter to use
%n = 100;
for i = 1:length(a)
    %i, a(i)
    if strcmp(a(i).Geometry,'Polygon')
        p_count = p_count + 1;
        
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
        %filter out the NaN
        latnan = find(isnan(a(i).Lat));
        a(i).Lat(latnan) = [];
        a(i).Lon(latnan) = [];
        lonnan = find(isnan(a(i).Lon));
        a(i).Lat(lonnan) = [];
        a(i).Lon(lonnan) = [];
        if length(a(i).Lon) ~= length(a(i).Lat)
            fprintf('Size mismatch between Lon and Lat\n')
        end
        if length(a(i).Lat) > n
            decimate = round(length(a(i).Lat)/n);
            lats = a(i).Lat(1:decimate:end);
            lons = a(i).Lon(1:decimate:end);
        else
            lats = a(i).Lat;
            lons = a(i).Lon;
        end
        
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
        
        %store perimter structure
        nums = length(lats);
        p_struct(p_count).det = [0 0 1 1 1];
        p_struct(p_count).power = ones(1,nums)*50;
        p_struct(p_count).data = ones(1,nums)*9;
        p_struct(p_count).conf = ones(1,nums)*95;
        p_struct(p_count).time = a(i).p_time;
        p_struct(p_count).lat = lats;
        p_struct(p_count).lon = lons;
        p_struct(p_count).file = replace(a(i).Name,' ','_');
        %p_struct(p_count).Name = a(i).Name;
        p_struct(p_count).xlon = [];
        p_struct(p_count).xlat = [];
        p_struct(p_count).fxdata = [];
        p_struct(p_count).axis = [];
        
    end
end %for
end % function
