function filled_countours(ps,an)
% plots filled contours for forecast and analysis

figure,contourf(ps.red.fxlong,ps.red.fxlat,ps.red.tign,20,'k')
hold on
scatter(ps.points(:,2),ps.points(:,1),'*r')
xlabel('Lon'),ylabel('Lat'),title('Forecast and active fire detections')
hold off

figure,contourf(ps.red.fxlong,ps.red.fxlat,an,20,'k')
hold on
scatter(ps.points(:,2),ps.points(:,1),'*r')
xlabel('Lon'),ylabel('Lat'),title('Analysis and active fire detections')
hold off

end