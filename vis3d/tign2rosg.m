function r=tign2ros(p,i)
[gx,gy]=grad(p.fxlong(:,:,i),p.fxlat(:,:,i),p.tign_g(:,:,i));
gn=sqrt(gx.^2+gy.^2);
r = 1./gn;
end