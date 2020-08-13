function ac = area_compare(w_a,w_b)
%compares areas of fires over time
% ts = '2013-08-17_07:00:00'
% w_a = read_wrfout_tign(wrf_a,ts);
% w_b = read_wrfout_tign(wrf_b,ts);

% w_a = read_wrfout_tign(wrf_a);
% w_b = read_wrfout_tign(wrf_b);

t_1 = max(min(w_a(:)),min(w_b(:)));
t_2 = min(max(w_a(:)),max(w_b(:)));

pts = 20;
t = linspace(t_1,t_2,pts);

for i = 1:pts
    area_a(i) = sum(sum(w_a < t(i)));
    area_b(i) = sum(sum(w_b < t(i)));
end

figure
plot(t/(24*3600),area_a,t/(24*3600),area_b)
legend('wrf-a','wrf-b')
%legend('Wet fuel','Normal Fuel')
xlabel('Time')
ylabel('Fire Area')
title('Comparison of Areas')
%xlim([0 6.1])
ac.area_a = area_a;
ac.area_b = area_b;

end



