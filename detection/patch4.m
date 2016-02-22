% to create conus.kml:
% download http://firemapper.sc.egov.usda.gov/data_viirs/kml/conus_hist/conus_20120914.kmz
% and gunzip 
% 
% to create w.mat:
% run Adam's simulation, currently results in
%
% /share_home/akochans/WRF341F/wrf-fire/WRFV3/test/em_utfire_1d_med_4km_200m
% then in Matlab
% arrays needed only once
% f='wrfout_d01_2013-08-20_00:00:00'; 
% t=nc2struct(f,{'Times'},{});  n=size(t.times,2)  
% w=nc2struct(f,{'Times','TIGN_G','FXLONG','FXLAT','UNIT_FXLAT','UNIT_FXLONG','XLONG','XLAT','NFUEL_CAT'},{},n);
% save ~/w.mat w    
%
% array at fire resolution every saved timestep
% to create s.mat:
% a=dir('wrfout_d01*');
% s=read_wrfout_sel({a.name},{'FGRNHFX',Times}); 
% save ~/s.mat s 
% 
% arrays at atm resolution every saved timestep
% to create ss.mat
% a=dir('wrfout_d01*')
% s=read_wrfout_sel({a.name},{'Times','UAH','VAH'})
% save ss s
% 
% fuels.m is created by WRF-SFIRE at the beginning of the run

% ****** REQUIRES Matlab 2013a - will not run in earlier versions *******

% run patch_load first

% figures
figmap=1;

% convert tign_g to datenum 
w.time=datenum(char(w.times)');
red.tign=(red.tign_g - max(red.tign_g(:)))/(24*60*60) + w.time;
min_tign=min(red.tign(:));
max_tign=max(red.tign(:));

cmap2=cmap;
cmap2(1:7,:)=NaN;
[cmap,imax]=cmapmod14;

figure(figmap);clf
lastdet=1;
for step=2:length(ss.time)  % over WRF frames
  figure(figmap);clf;hold off
  for ipass=1:2,    
    det=find(r.time <= ss.time(step));
    ndet=length(det);
    % detections to now
    for idet=1:ndet,
        x=r.x{det(idet)}; % load fire detection image 
        age=t-r.time(idet); % age of detection in days
        offset = min(imax,floor(4*age)); % offset in colormap for detection age
        dd=x.data(:)>6;  % indices of detected
        x.data(dd)=x.data(dd)+3*offset;   % transition to yellow
        fprintf('step %i pass %i granule %i detections %i\n',step,ipass,idet,sum(dd))
%        if idet<ndet, % all but the last
%            if any(dd),
%                fprintf(' %i fire detections, rest transparent\n',nnz(dd));
%                d0=find(x.data(:)<=6);
%                x.data(d0)=0;
%                showmod14(x)
%                hold on
%            else
%                fprintf(' all transparent, skipping\n')
%            end
%        else
%            fprintf(' %i fire detections, displaying granule\n',nnz(dd))
%            showmod14(x)
%            hold on
%        end
        if ipass==1, % build background
               showmod14(x)
        elseif any(dd), % all except fire transparent
               x.data(~dd)=0;
               showmod14(x)
        end
        hold on
      end
    end
    u=ss.uh(:,:,step);
    v=ss.vh(:,:,step);
    maxw=max(sqrt(u(:).*u(:)+v(:).*v(:)));
    fprintf('step %i max windspeed %g granules %i\n',step,maxw,ndet)
    sc=0.006;quiver(w.xlong,w.xlat,sc*u,sc*v,0); % wind
    hold on
    t=ss.time(step);
    contour(red.fxlong,red.fxlat,red.tign,[t t],'-k'); % fireline
    title(datestr(t));   
    axis(red.axis)
    hold off
    drawnow
    pause(0.1)
    M(step)=getframe(gcf);
    
    title(datestr(t));
end

mov2avi(M,'Mw')