function r=lfn2ros(p,i)
lfn=p.lfn(:,:,i);
lfn_start=p.lfn_start(:,:,i);
fxlong=p.fxlong(:,:,i);
fxlat=p.fxlat(:,:,i);
dt=p.dt;
[gx,gy]=grad(fxlong,fxlat,lfn);
gn=sqrt(gx.^2+gy.^2);
dLdt = (lfn_start-lfn)/dt;
r=dLdt(2:end-1,2:end-1)./gn;
end