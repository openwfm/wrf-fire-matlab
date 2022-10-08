function r=tign2ros(p,i)
dxlon=(p.fxlong(3:end,2:end-1,i)-p.fxlong(1:end-2,2:end-1,i));
dxlat=(p.fxlat(3:end,2:end-1,i)-p.fxlong(1:end-2,2:end-1,i));
dx=sqrt(dxlon.^2+dxlat.^2);
gx=(p.tign_g(3:end,2:end-1,i)-p.tign_g(1:end-2,2:end-1,i))./dx;
dylon=(p.fxlong(2:end-1,3:end,i)-p.fxlong(2:end-1,1:end-2,i));
dylat=(p.fxlat(2:end-1,3:end,i)-p.fxlong(2:end-1,1:end-2,i));
dy=sqrt(dylon.^2+dylat.^2);
gy=(p.tign_g(2:end-1,3:end,i)-p.tign_g(2:end-1,1:end-2,i))./dy;
gn=sqrt(gx.^2+gy.^2);
r = 1./gn;
end