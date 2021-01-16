function filled_countours(ps,an)
% plots filled contours for forecast and analysis times in seconds

tmax = max(ps.red.tign_g(:))-3600;
tmin = min(ps.red.tign_g(:));
clines = linspace(tmin,tmax,20);

if isfield(ps.red,'red')
    figure,contourf(ps.red.red.fxlong(1:10:end,1:10:end),ps.red.red.fxlat(1:10:end,1:10:end) ...
        ,ps.red.red.tign(1:10:end,1:10:end),clines,'k')
else
    figure,contourf(ps.red.fxlong,ps.red.fxlat,ps.red.tign,clines,'k')
end
hold on
scatter(ps.points(:,2),ps.points(:,1),'*r')
xlabel('Lon'),ylabel('Lat'),title('Forecast and active fire detections')
hold off


figure,contourf(ps.red.fxlong,ps.red.fxlat,an,clines,'k')
hold on
scatter(ps.points(:,2),ps.points(:,1),'*r')
xlabel('Lon'),ylabel('Lat'),title('Analysis and active fire detections')
hold off

end