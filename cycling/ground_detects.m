function gd = ground_detects(red)

if exist('g_ground.mat','file')
    load g_ground.mat
else
    gs = l2_detect(red);
    save g_ground.mat gs
end
%fires = gs.fires;
g = gs.g;
n = length(g);
water = [];
land = [];
for i = 1:n
    wat_msk = g(i).data == 3;
    lan_msk = g(i).data == 5;
    wat_time = g(i).time*ones(sum(wat_msk),1);
    lan_time = g(i).time*ones(sum(lan_msk),1);
    water = [water;[g(i).lat(wat_msk),g(i).lon(wat_msk),wat_time]];
    land= [land;[g(i).lat(lan_msk),g(i).lon(lan_msk),lan_time]];
end
land = [fixpoints2grid(red,land),land(:,3)];
water = [fixpoints2grid(red,water),water(:,3)];
land = unique(land,'rows');
water = unique(water,'rows');
gd.water = water;
gd.land = land;


end
