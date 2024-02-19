function r=tign2rosg(p,i)
t=p.tign_g(:,:,i);
t(t==t(1,1))=nan;
[gx,gy]=grad(p.fxlong(:,:,i),p.fxlat(:,:,i),t);
gn=sqrt(gx.^2+gy.^2);
r = 1./gn;
end