function cone_compare(ps,tign2)
%compares feature of fire arrival cones
%ps = tign_try(w), tign = squish(ps)



lon = ps.red.fxlong;
lat = ps.red.fxlat;
tign = ps.red.tign;
%blur the data for smoother gradients
tign = imgaussfilt(tign,1);
tign2 = imgaussfilt(tign2,1);

%compute are of fires
t_end = min(max(tign(:)),max(tign2(:)));
a1 = sum(sum(tign<t_end));
a2 = sum(sum(tign2<t_end));
fprintf('Forecast Area: %d Data area: %d \n',a1,a2);

cull = input_num('Thin data set? [1]',1);
lon = ps.red.fxlong(1:cull:end,1:cull:end);
lat = ps.red.fxlat(1:cull:end,1:cull:end);
tign = ps.red.tign(1:cull:end,1:cull:end);
tign2 = tign2(1:cull:end,1:cull:end);

%compute gradient step sizes
E = wgs84Ellipsoid;
[n,m] = size(lon);
hx = distance(lat(1,round(m/2)),lon(1,1),lat(1,round(m/2)),lon(end,end),E)/m;
hy = distance(lat(1,1),lon(round(n/2),1),lat(end,end),lon(round(n/2),1),E)/n;

[dx1,dy1] = fire_gradients(lon,lat,tign,1);
dx1=dx1/hx;dy1=dy1/hy;
[dx2,dy2] = fire_gradients(lon,lat,tign2,1);
dx2=dx2/hx;dy2=dy2/hy;

%mask for only the fire cone
t_msk1 = tign<max(tign(:));
t_msk2 = tign2<max(tign2(:));
t_msk = logical(t_msk1.*t_msk2);

%measure angles
theta1 = atan2(dy1(t_msk),dx1(t_msk));
theta2 = atan2(dy2(t_msk),dx2(t_msk));

td = theta1-theta2;
td(td>pi) = td(td>pi)-2*pi;
td(td<-pi) = td(td<-pi)+2*pi;
td = td(~isnan(td));
figure,histogram(td)
format short
tstr= sprintf('Angle difference in gradients \n Mean : %f Std deviation: %f',mean(td),std(td));
title(tstr)
xlabel('Difference of angles (radians)')
ylabel('Number')

figure
quiver(lon(t_msk),lat(t_msk),dy1(t_msk),dx1(t_msk))
% quiver(lon,lat,dy1,dx1)
hold on
quiver(lon(t_msk),lat(t_msk),dy2(t_msk),dx2(t_msk))
% quiver(lon,lat,dy2,dx2)
title('Gradients in fire surfaces, unit vectors')
legend('Forecast','Interpolated')

%get vectors which are not unit vectors
[du1,dv1] = fire_gradients(lon,lat,tign,0);
[du2,dv2] = fire_gradients(lon,lat,tign2,0);
du1=du1/hx;dv1=dv1/hy;
du2=du2/hx;dv2=dv2/hy;


figure
quiver(lon(t_msk),lat(t_msk),dv1(t_msk),du1(t_msk))
hold on
quiver(lon(t_msk),lat(t_msk),dv2(t_msk),du2(t_msk))
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


%rate of spread
rx1 = 1./du1/(24*3600);
ry1 = 1./dv1/(24*3600);
rx2 = 1./du2/(24*3600);
ry2 = 1./dv2/(24*3600);
r1 = sqrt(rx1.^2+ry1.^2);
r2 = sqrt(rx2.^2+ry2.^2);

r_diff = r1-r2;
r_msk = abs(r_diff)<10;
%r_diff = r_diff(abs(r_diff)<10);
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

figure,histogram(ps.red.nfuel_cat(r_fast)),xticks(1:14)
title('Fuel Types Where Fire is Burning too Fast')
xlabel('Fuel Type'),ylabel('Number')
figure,histogram(ps.red.nfuel_cat(r_slow)),xticks(1:14)
title('Fuel Types Where Fire is Burning too Slow')
xlabel('Fuel Type'),ylabel('Number')
% figurelow
% quiver(lon(t_msk),lat(t_msk),rx1(t_msk),ry1(t_msk))
% hold on
% quiver(lon(t_msk),lat(t_msk),rx2(t_msk),ry2(t_msk))
% title('ROS vectors')
% legend('Ground Truth','Interpolated')

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



end