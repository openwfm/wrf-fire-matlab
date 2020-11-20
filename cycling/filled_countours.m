function filled_countours(ps,an)
% plots filled contours for forecast and analysis

if isfield(ps.red,'ornge')
    figure,contourf(ps.red.red.fxlong(1:10:end,1:10:end),ps.red.red.fxlat(1:10:end,1:10:end) ...
        ,ps.red.red.tign(1:10:end,1:10:end),20,'k')
else
    figure,contourf(ps.red.fxlong,ps.red.fxlat,ps.red.tign,20,'k')
end
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