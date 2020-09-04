function ac = area_compare(ps,tign_b)
%compares areas of fires over time
% outputs graph, areas and a mask
% ps = cluster_paths(w) --- Forecast
% tign_b = quish2(ps) --- "Analysis"

tign_a = ps.red.tign;
%create a mask 
msk=tign_b-tign_a;
msk(abs(msk)<0.2) = 0;
msk(msk>0)=1; %add fuel moisture ==> slow fire
msk(msk<0)=-1; %subtract fuel moist ==> speed up fire

%time +- 1/4 day at simulation start and end
t_1 = max(min(tign_a(:)),min(tign_b(:)))+0.25;
t_2 = min(max(tign_a(:)),max(tign_b(:)))-0.25;

pts = 20;
t = linspace(t_1,t_2,pts);

for i = 1:pts
    area_a(i) = sum(sum(tign_a < t(i)));
    area_b(i) = sum(sum(tign_b < t(i)));
end

figure
plot(t/(24*3600),area_a,t/(24*3600),area_b)
%legend('wrf-a','wrf-b')
legend('Forecast','Estimate')
%legend('Wet fuel','Normal Fuel')
xlabel('Time')
ylabel('Fire Area')
title('Comparison of Areas')

% figure,mesh(ps.red.fxlong,ps.red.fxlat,msk)
% title('Mask for fuel adjustments')
% xlabel('lon'),ylabel('lat')
ac.area_a = area_a;
ac.area_b = area_b;
ac.msk = msk;

end



