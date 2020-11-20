function plot_script(red)
% script plots all detections 

%red = red_ps;
hold on
tign_days = red.tign_g/(24*3600);


%plot_state(1,red,'Forecast from cycling',red.tign_g,g,time_bounds(1:2));
%hold on
figure,contourf(red.fxlong,red.fxlat,tign_days,20);
%figure,mesh(red.fxlong,red.fxlat,tign_days);
xlabel('Lon'),ylabel('Lat')
%zlabel('Time [days]')
%hold off
end