function plot_paths2d(path_struct)
%plots paths in 2d setting
% input: path_struct - structure with graph, paths, red, etc...
%        path_struct = graph_dets(w,cull)
r = path_struct.red;

% draw filled contour of fuels
fuel_map = r.nfuel_cat;
%fuel_map = simp_fuels(fuel_map);
figure,contourf(r.fxlong,r.fxlat,fuel_map,'LineWidth',0.5);
colormap summer;

%draw paths on figure
paths = path_struct.paths;
pts = path_struct.points;
hold on
for i = 1:length(paths)
    p = paths(i).p;
    plot(pts(p,2),pts(p,1),'r','LineWidth',2);
end %for i
hold off

%plot elevation data and paths

% draw filled contour of fuels
elev = r.fhgt;
%fuel_map = simp_fuels(fuel_map);
figure,contour(r.fxlong,r.fxlat,elev,'k','LineWidth',0.5,'ShowText','on');

%draw paths on figure
paths = path_struct.paths;
pts = path_struct.points;
hold on
for i = 1:length(paths)
    p = paths(i).p;
    plot(pts(p,2),pts(p,1),'r','LineWidth',2);
end %for i
hold off

%plot paths with gradient of surface
%smooth surface first
r.tign = imgaussfilt(r.tign,1/2)
[dx,dy]= fire_gradients(r.fxlong,r.fxlat,r.tign,1);

figure,contour(r.fxlong,r.fxlat,r.tign)
hold on
quiver(r.fxlong,r.fxlat,dy,dx)
for i = 1:length(paths)
    p = paths(i).p;
    plot(pts(p,2),pts(p,1),'r','LineWidth',1);
end %for i
hold off


end %function