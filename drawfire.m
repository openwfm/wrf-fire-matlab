function p=drawfire(f)
% p=drawfire(f)
% show heatflux stored in wrfout
vars={'FIRE_AREA','FGRNHFX','LFN','TIGN_G','FUEL_FRAC','FXLONG','FXLAT','NFUEL_CAT','Times'};
p=nc2struct(f,vars,{'DX','DY'});
nframes=size(p.fgrnhfx,3);
[i,j,a]=find(p.fgrnhfx(:,:,end));
is=min(i);
ie=max(i);
js=min(j);
je=max(j);
r=max(ie-is,je-js)/2;
ci = (is+ie)/2;
cj = (js+je)/2;
is=round(ci-r);
ie=round(ci+r);
js=round(cj-r);
je=round(cj+r);
for n=2:nframes
    figure(1)
    subplot(1,2,1)
    mesh(p.fxlong(is:ie,js:je,n),p.fxlat(is:ie,js:je,n),p.fgrnhfx(is:ie,js:je,n)),
    title(['FGRNHFX frame ',num2str(n)])
    subplot(1,2,2)
    %{
    mesh(p.fxlong(is:ie,js:je,n),p.fxlat(is:ie,js:je,n),p.fire_area(is:ie,js:je,n))
    hold on
    contour(p.fxlong(is:ie,js:je,n),p.fxlat(is:ie,js:je,n),p.lfn(is:ie,js:je,n),[0 0],'k')
    hold off
    view(0.9,90)
    title(['FIRE\_AREA frame ',num2str(n)])
    subplot(2,2,3)
    mesh(p.fxlong(is:ie,js:je,n),p.fxlat(is:ie,js:je,n),p.fire_area(is:ie,js:je,n-1))
    view(0.9,90)
    title(['FIRE\_AREA frame ',num2str(n-1)])
    subplot(2,2,4)
    %}
    mesh(p.fxlong(is:ie,js:je,n),p.fxlat(is:ie,js:je,n),...
        p.fuel_frac(is:ie,js:je,n-1)-p.fuel_frac(is:ie,js:je,n))
    title(['FUEL\_FRAC frame ',num2str(n),' decrement'])
    drawnow
    pause(0.3)
    % figure(2)
    % mesh(p.fxlong(is:ie,js:je,n),p.fxlat(is:ie,js:je,n),p.fire_area(is:ie,js:je,n))
    % view(0.9,90)
end
