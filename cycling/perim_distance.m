function [d,v] =perim_distance(w,g)
%computes distance from ignition to furthest point on perimeter, avg ROS to get there 
%in straight line
% input   -   w, struct with tign
%                  g, struct with detections
 %output   -d,v     distance (m), ros (m/s)
 
 r = subset_domain(w);
 min_t = min(r.tign_g(:));
 max_t = max(r.tign_g(:)-10);
 
 t = max_t-min_t;
 [x0,y0] = find(r.tign_g==min_t);
 fm = r.tign_g < max_t;
 longs = r.fxlong(fm);
 longs = longs(:);
 lats = r.fxlat(fm);
 lats = lats(:);
 
 E = wgs84Ellipsoid;
 %loop through the distances...
 point = [r.fxlat(x0,y0),r.fxlong(x0,y0)];
 ignition_point = point;
 d = 0;
 
 for i=1:10:length(longs)
     new_point = [lats(i),longs(i)];
     new_d = distance(new_point,ignition_point,E);
     if new_d > d
         d = new_d;
         point = new_point
     end
 end
 v = d/t;
 
 dd = 0;
 det_ignition = [0 0];
 %compute distance in satellite data
 for j = 1:length(g)
     if sum(g(j).det(3:5)) > 0
         dm = g(j).data > 6;
         det_longs = g(j).xlon(dm);
         det_longs = det_longs(:);
         det_lats = g(j).xlat(dm);
         det_lats = det_lats(:);
         %pick ignition point
         if norm(det_ignition) ==0
             det_ignition(1) = mean(det_lats);
             det_ignition(2) = mean(det_longs);
             det_start = g(j).time;
         else 
             for k = 1:length(det_lats)
                 new_det_point = [det_lats(k),det_longs(k)];
                 new_dd = distance(new_det_point,det_ignition,E);
                 if new_dd > dd
                     dd = new_dd;
                     det_point = new_det_point;
                     det_time = g(j).time;
                 end
                 
             end
         end
     end
 end
 det_t = (det_time-det_start)*(3600*24);
 dv = dd/det_t;
 fprintf('Forecast distanc: %0.2f     Forecast ROS: %0.2f  \n',d,v)
 fprintf('Detection distance : %0.2f   Detection ROS: %0.2f \n',dd,dv);
 
 n = 10;
 c_lines = linspace(min_t,max_t,n);
 
 figure,contour(r.fxlong,r.fxlat,r.tign_g,c_lines,'k')
 hold on
 scatter(ignition_point(2),ignition_point(1),500,'*r')
 scatter(point(2),point(1),500,'r*');
 title('Foecast points')
 hold off
 
 figure,contour(r.fxlong,r.fxlat,r.tign_g,c_lines,'k')
 title('Detection Points');
 hold on
 scatter(det_ignition(2),det_ignition(1),500,'*r');
 scatter(det_point(2),det_point(1),500,'*r');
 hold off
 
end