function perim_plot(t1,t2,pc)
%plots final perimeter of two simulations
%inputs t1, t2 - tign in matrix, 
%        pc - path structure pc = cluster_paths(w,1)

mxt = max(t1(:))-0.2;
p = [mxt mxt];
%p = 20;

figure,contour(pc.red.fxlong,pc.red.fxlat,t1,p,'k');
hold on
contour(pc.red.fxlong,pc.red.fxlat,t2,p,'r')
title('Perimeters')
legend('t1','t2')


end % function
