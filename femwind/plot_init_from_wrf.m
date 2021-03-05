p=nc2struct(f,{'HT_FMW','U0_FMW','V0_FMW','W0_FMW','FXLONG','FXLAT','ZSF',...
    'UNIT_FXLAT','UNIT_FXLONG','XLONG','XLAT'},{'DX','DY'},60)
[m,n]=size(p.xlat);
[fm,fn]=size(p.fxlat);
fdx=p.dx*(m/fm);
fdy=p.dy*(n/fn);
n=size(p.u0_fmw);
CX{1}=zeros(n);
CX{2}=zeros(n);
CX{3}=zeros(n);
x = fdx*([1:n(1)]-0.5);
y = fdy*([1:n(2)]-0.5);
[xx,yy]=ndgrid(x,y);
for k=1:n(3)
    CX{1}(:,:,k)=xx;
    CX{2}(:,:,k)=yy;
    CX{3}(:,:,k)=p.zsf+p.ht_fmw(k);
end
W={p.u0_fmw,p.v0_fmw,p.w0_fmw};
plot_wind_3d(CX,W,[],[],20)