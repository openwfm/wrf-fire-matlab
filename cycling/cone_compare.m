function cone_compare(ps,tign2)
%compares feature of fire arrival cones
%ps = tign_try(w), tign = squish(ps)

lon = ps.red.fxlong;
lat = ps.red.fxlat;
%forecast
tign = ps.red.tign;
%blur the data for smoother gradients
% st = 3;
% tign = imgaussfilt(tign,st);
% tign2 = imgaussfilt(tign2,st);

%compare areas of the two and plot over time
area_compare(ps,tign2); 

%compute area of fires
t_end = min(max(tign(:)),max(tign2(:)))-0.1;
a1 = sum(sum(tign<t_end));
a2 = sum(sum(tign2<t_end));
%make the top flate for each
% tign(tign>=t_end)=t_end;
% tign2(tign2>=t_end)=t_end;
fprintf('Forecast Area: %d Data area: %d \n',a1,a2);
figure,contour(lon,lat,tign,[t_end t_end],'k')
hold on,contour(lon,lat,tign2,[t_end t_end],'b')
legend('forecast','estimate')
t_str = sprintf('Perimeters \n Forecast Area = %d    Data Area = %d',a1,a2);
title(t_str)
cull = input_num('Thin data set? [1]',1,1);
if cull ~= 1
lon = lon(1:cull:end,1:cull:end);
lat = lat(1:cull:end,1:cull:end);
tign = tign(1:cull:end,1:cull:end);
tign2 = tign2(1:cull:end,1:cull:end);
end

%compute gradient step sizes
E = wgs84Ellipsoid;
[n,m] = size(lon);
hx = distance(lat(1,round(m/2)),lon(1,1),lat(1,round(m/2)),lon(end,end),E)/m;
hy = distance(lat(1,1),lon(round(n/2),1),lat(end,end),lon(round(n/2),1),E)/n;
aspect_ratio = [1 hy/hx 1];


%gradient in fire arrival time, unit vector
[dx1,dy1] = fire_gradients(lon,lat,tign,1);
% [dx1,dy1] = gradient(tign);
% dx = dx1./sqrt(dx1.^2+dy1.^2);
% dy = dy1./sqrt(dx1.^2+dy1.^2);
% dx1=dx;
% dy1=dy;
% clear dx dy
figure,contour(lon,lat,tign,20,'k');
daspect(aspect_ratio)
hold on
quiver(lon,lat,dx1,dy1)
hold off
% dx1=dx1/hx;dy1=dy1/hy;
%unit vector
[dx2,dy2] = fire_gradients(lon,lat,tign2,1);
% figure,contour(lon,lat,tign2,20,'k');
% hold on
% quiver(lon,lat,dx2,dy2)
% hold off
%[dx2,dy2] = gradient(tign2);
% dx2=dx2/hx;dy2=dy2/hy;

%compute the gradient of the terrain
elev = ps.red.fhgt;
[aspect,slope,ey,ex] = gradientm(lat,lon,elev,E);
% figure,contour(lon,lat,elev,'k')
% hold on
% quiver(lon(1:5:end,1:5:end),lat(1:5:end,1:5:end),ex(1:5:end,1:5:end),ey(1:5:end,1:5:end))

%slopes of fire directions by directional derivatives
sl1 = ex.*dx1+ey.*dy1;
sl2 = ex.*dx2+ey.*dy2;
sl_diff = sl1-sl2;
figure,histogram(sl_diff)
t_str=sprintf('Difference in slopes normal to fire front \n Mean = %f', mean(sl_diff(~isnan(sl_diff))));
xlabel('Slope differences')
ylabel('Number')
title(t_str)

%mask for only the fire cone
t_msk1 = tign<max(tign(:))-0.1;
t_msk2 = tign2<max(tign2(:))-0.1;
t_msk = logical(t_msk1.*t_msk2);

%measure angles
% theta1 = atan2(dy1(t_msk),dx1(t_msk));
% theta2 = atan2(dy2(t_msk),dx2(t_msk));
%no mask
theta1 = atan2(dy1,dx1);
theta2 = atan2(dy2,dx2);

td = theta1-theta2;
td(td>pi) = td(td>pi)-2*pi;
td(td<-pi) = td(td<-pi)+2*pi;
td_msk = ~isnan(td);
b_msk = abs(td)<pi/6;
figure,histogram(td)
format short
tstr= sprintf('Angle difference in gradients \n Mean : %f Std deviation: %f',mean(td(td_msk)),std(td(td_msk)));
title(tstr)
xlabel('Difference of angles (radians)')
ylabel('Number')

figure
quiver(lon(t_msk),lat(t_msk),dx1(t_msk),dy1(t_msk))
% quiver(lon,lat,dy1,dx1)
hold on
quiver(lon(t_msk),lat(t_msk),dx2(t_msk),dy2(t_msk))
% quiver(lon,lat,dx2,dy2)
title('Gradients in fire surfaces, unit vectors')
legend('Forecast','Interpolated')

%get vectors which are not unit vectors
[du1,dv1] = fire_gradients(lon,lat,tign,0);
[du2,dv2] = fire_gradients(lon,lat,tign2,0);

figure
quiver(lon(t_msk),lat(t_msk),du1(t_msk),du1(t_msk))
hold on
quiver(lon(t_msk),lat(t_msk),du2(t_msk),dv2(t_msk))
title('Gradients in fire surfaces')
legend('Forecast','Interpolated')

%get magnitudes of these vectors for comparison
% m1 = sqrt(du1(t_msk).^2+dv1(t_msk).^2);
% m2 = sqrt(du2(t_msk).^2+dv2(t_msk).^2);
m1 = sqrt(du1.^2+dv1.^2);
m2 = sqrt(du2.^2+dv2.^2);
mdiff = m1-m2;
g_msk = abs(mdiff)>=0.0;
% figure,histogram(mdiff(g_msk));
% avg_mdiff = mean(mdiff(g_msk));
% std_mdiff = std(mdiff(g_msk));
% tstr = sprintf('Histogram of difference in vector magnitudes \n mean = %f std = %f',avg_mdiff,std_mdiff);
% title(tstr);

%use the slope from gradientm function instead
%rate of spread
rx1 = 1./du1/(24*3600);
ry1 = 1./dv1/(24*3600);
rx2 = 1./du2/(24*3600);
ry2 = 1./dv2/(24*3600);
r1 = sqrt(rx1.^2+ry1.^2);
r2 = sqrt(rx2.^2+ry2.^2);
cut_off = est_max(ps,r2)
%cut_off = 0.1;
b1 = r1<cut_off;b2 = r2 < cut_off;b_msk = logical(b1.*b2);
figure
quiver(lon(b_msk),lat(b_msk),rx1(b_msk),ry1(b_msk))
hold on
quiver(lon(b_msk),lat(b_msk),rx2(b_msk),ry2(b_msk))
t_str =sprintf('ROS vectors for ROS < %f',cut_off);
title(t_str);

r_diff = r1-r2;
%r_diff = imgaussfilt(r_diff,1/2);
r_msk = abs(r_diff)<10;
avg_r_diff = mean(r_diff(r_msk));
std_r_diff = std(r_diff(r_msk));
figure,histogram(r_diff(r_msk));
tstr= sprintf('Differences in ROS, forecast-estimate \n Mean = %f  Std Dev. = %f',avg_r_diff,std_r_diff);
title(tstr)
xlabel('ROS (m/s)')
ylabel('Number')
r_fast = r_diff>0;%(avg_r_diff+1*std_r_diff);
r_slow = r_diff<0;%(avg_r_diff-1*std_r_diff);
figure,scatter(lon(r_fast),lat(r_fast),'*r');
hold on,scatter(lon(r_slow),lat(r_slow),'b');
title('Locations for fuel adjustment')
legend('Forecast too fast','Forecast too slow')

%regression on slope and ros differences
r_diff(abs(r_diff)>2) = NaN;
sl_diff(abs(r_diff)>2) = NaN;
figure,scatter(sl_diff(~isnan(r_diff)),r_diff(~isnan(r_diff)));
mdl = fitlm(sl_diff(:),r_diff(:));

% figure,histogram(ps.red.nfuel_cat(r_fast)),xticks(1:14)
% title('Fuel Types Where Fire is Burning too Fast')
% xlabel('Fuel Type'),ylabel('Number')
% figure,histogram(ps.red.nfuel_cat(r_slow)),xticks(1:14)
% title('Fuel Types Where Fire is Burning too Slow')
% xlabel('Fuel Type'),ylabel('Number')

% figure
% quiver(lon(t_msk),lat(t_msk),rx1(t_msk),ry1(t_msk))
% hold on
% quiver(lon(t_msk),lat(t_msk),rx2(t_msk),ry2(t_msk))
% title('ROS vectors')
% legend('Ground Truth','Interpolated')
figure,histogram(ps.red.nfuel_cat(t_msk)),xticks(1:14)
title('Fuel types')


%find where the fire is burning too fast

% fast = mdiff < avg_mdiff-1/2*std_mdiff;
% slow = mdiff > avg_mdiff+1/2*std_mdiff;
%figure,scatter(lon(slow),lat(slow)),
% figure,scatter(lon(fast),lat(fast),'*r')
%legend('Forecast is Slow','Forecat is fast')

%cluster the domain by fuel type and gradient difference
% [s_idx,s_c] = kmeans([lon(:),lat(:),mdiff(:)],2);
% figure,scatter(lon(s_idx==1),lat(s_idx==1),'r')
% hold on,scatter(lon(s_idx==2),lat(s_idx==2),'b')
% legend('Forecast too fast','Forecast too slow')

%compute slope of terrain

%compare by fuel type
for i = 1:13
    fuel_mask = ps.red.nfuel_cat == i;
    msk = logical(fuel_mask.*r_msk);
    fuel_rate_diff =  r1(msk)-r2(msk);
    mean_frd = mean(fuel_rate_diff);
    std_frd = std(fuel_rate_diff);
    format short
    fprintf('Fuel type: %d ROS difference: %f Std. Deviation: %f \n',i,mean_frd,std_frd);
    
end




end