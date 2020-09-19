function sm_up = smooth_up(lon,lat,tign)

t0 = min(tign(:));
t1 = max(tign(:));

steps = 20;
t = linspace(t0,t1,steps);

%compute different std. for imgausfilt
mx = 2;
mn = 1/2;
m = (mn-mx)/(steps-2);
y = @(t) mx + m*(t-2);

%keep track of the change in ignition time
t_diff(1) = 0;

t_temp = tign;
% figure(136),mesh(lon,lat,tign),title('Original')
% figure(137),mesh(lon,lat,tign),title('Smoothing')
for i = 2:steps
    m1 = tign < t(i-1);
    m2 = tign < t(i);
    msk = logical(m2-m1);
    t_blur = imgaussfilt(tign,y(i));
    t_temp(msk) = t_blur(msk);
    figure(137),mesh(lon,lat,t_temp)
    pause(5/steps)
    t_diff(i) = min(t_temp(:))-t0;    
    t_temp(m2) = t_temp(m2)-i/steps*t_diff(i);
    %t_temp = t_temp-i/steps*t_diff(i);
end
    %shift down
%     pixel_time_diff = t_temp-tign;
%     time_slope = pixel_time_diff./t_diff(i);
%     t_temp3 = t_temp - time_slope.*pixel_time_diff;
%     t_shift = max(t_temp3,t_temp);
% figure,plot(t_diff),title('Change in tign')
sm_up = t_temp;