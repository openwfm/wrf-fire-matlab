function w=read_wrfout_tign(f,ts)
    % w=read_wrfout_tign(f) 
    % w=read_wrfout_tign(f,timestep) 
    % read variables from wrfout f for detect_fit_level2 and cycles
    % in
    %    f file name
    %    ts time step (string, optional). Read last time step if not given.
    % 
    
    if f(end) == 't'
        fprintf('Matlab file as input, loading \n')
        load(f)
        return
    end
    t=nc2struct(f,{'Times'},{});  nframes=size(t.times,2);
    alltimes=char(t.times')
    fprintf('Last time step in %s is %i at %s\n',f,nframes,alltimes(nframes,:))
    if exist('ts','var')
        frame=0;
        for i=nframes:-1:1,
            if strcmp(ts,alltimes(i,:))
                fprintf('Found timestep %i at %s\n',i,ts)
                frame=i;
                break

            end
        end
%         if ~isstring(ts)
%             frame = ts;
%         else
%             frame = input_num('Enter frame = ? ',nframes)
%         end
        fprintf('Choosing frame %d \n',frame)
        if(frame==0),
            warning(['Time step ',ts,' not found'])
            w=[];
            return
        end
    else
        frame=nframes;
    end

    w=nc2struct(f,{'Times','FGRNHFX','TIGN_G','FXLONG','FXLAT','UNIT_FXLAT','UNIT_FXLONG','LFN','UF','VF',...
        'XLONG','XLAT','NFUEL_CAT','ITIMESTEP','FMC_G','ROS','HGT','FUEL_FRAC','FUEL_FRAC_BURNT'},{'DX','DY','DT'},frame);
    F = scatteredInterpolant(w.xlong(:),w.xlat(:),w.hgt(:));
    w.fhgt = F(w.fxlong,w.fxlat);
    %w.fhgt = griddata(w.xlong,w.xlat,w.hgt,w.fxlong,w.fxlat);
    w.times=char(w.times');
    w.cwd=pwd;
    w.datestr=datestr(clock);
    s=ls('-l',f);w.ls=s(1:end-1);
    s=ls('-lL',f);w.lsL=s(1:end-1);
end
