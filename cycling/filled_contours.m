function filled_countours(ps,an)
% plots filled contours for forecast and analysis times in seconds



tmax = max(ps.points(:,3));
tmax = ps.red.end_datenum;
tmin = min(min(ps.red.tign_g(:)),min(an(:)));
clines = linspace(tmin,tmax,20);

if isfield(ps.red,'red')
    figure,contourf(ps.red.red.fxlong(1:10:end,1:10:end),ps.red.red.fxlat(1:10:end,1:10:end) ...
        ,ps.red.red.tign(1:10:end,1:10:end),clines,'k')
else
    %figure,contourf(ps.red.fxlong,ps.red.fxlat,ps.red.tign_g,clines,'k')
    figure,contourf(ps.red.fxlong,ps.red.fxlat,ps.red.tign_g,'k')
end
hold on

points_set = ps.points(:,3)<=tmax;
scatter(ps.points(points_set,2),ps.points(points_set,1),'*r')
xlabel('Lon'),ylabel('Lat')
%title('Forecast and active fire detections')
title('"Ground Truth" and artificial fire detections')
hold off


%figure,contourf(ps.red.fxlong,ps.red.fxlat,an,clines,'k')
figure,contourf(ps.red.fxlong,ps.red.fxlat,an,'k')
hold on
scatter(ps.points(:,2),ps.points(:,1),'*r')
xlabel('Lon'),ylabel('Lat')
%title('Analysis and active fire detections')
title('Estimate and artificial fire detections')
hold off

end