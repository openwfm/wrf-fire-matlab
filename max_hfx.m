% function max_hfx(p,n)
% p from drawfire, n=frame number
[i,j,v]=find(p.fgrnhfx(:,:,n));
[ni,nj]=size(p.fgrnhfx(:,:,n));
[vm,km]=max(v);
km=km(1);
im=i(km);
jm=j(km);
fprintf('max fgrnhfx %f at (%i %i) = (%f %f) fxlong/fxlat\n',...
    p.fgrnhfx(im,jm,n),im,jm,p.fxlong(im,jm,n),p.fxlat(im,jm,n))
r=5;
is=max(1,im-r);
ie=min(ni,im+r);
js=max(1,jm-r);
je=min(nj,jm+r);
fprintf('displaying (%i:%i,%i:%i) of (1:%i,1:%i)\n',is,ie,js,ie,ni,nj)
format short
fgrnhfx=p.fgrnhfx(is:ie,js:je,n)
fuel_frac=p.fuel_frac(is:ie,js:je,n)
lfn=p.lfn(is:ie,js:je,n)
tign_g=p.tign_g(is:ie,js:je,n)
