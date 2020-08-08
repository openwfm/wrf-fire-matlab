function cone_compare(ps,tign2)
%compares feature of fire arrival cones
%ps = tign_try(w), tign = squish(ps)

lon = ps.red.fxlong;
lat = ps.red.fxlat;
tign = ps.red.tign;

[dx1,dy1] = fire_gradients(lon,lat,tign,1);
[dx2,dy2] = fire_gradients(lon,lat,tign2,1);

%mask for only the fire cone
t_msk = tign<max(tign(:));

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
hold on
quiver(lon(t_msk),lat(t_msk),dy2(t_msk),dx2(t_msk))
title('Gradients in fire surfaces, unit vectors')
legend('Ground Truth','Interpolated')

%get vectors which are not unit vectors
[du1,dv1] = fire_gradients(lon,lat,tign,0);
[du2,dv2] = fire_gradients(lon,lat,tign2,0);
figure
quiver(lon(t_msk),lat(t_msk),dv1(t_msk),du1(t_msk))
hold on
quiver(lon(t_msk),lat(t_msk),dv2(t_msk),du2(t_msk))
title('Gradients in fire surfaces')
legend('Ground Truth','Interpolated')




end