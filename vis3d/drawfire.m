function p=drawfire(f,ifire,nframes)
% p=drawfire(f)
% show heatflux stored in wrfout
if ischar(f)
    vars={'FIRE_AREA','FGRNHFX','LFN','TIGN_G','FUEL_FRAC','FXLONG','FXLAT','NFUEL_CAT','Times'};
    p=nc2struct(f,vars,{'DX','DY'});
    % all preprocessing here 
    p.times=char(p.times');
    for n=1:size(p.fgrnhfx,3)
        p.fgrnhfx_max(n)=max(max(p.fgrnhfx(:,:,n)));
    end
else
    p=f;  % structure passed in
end
% graphics only
if ~exist('nframes','var')
    nframes=size(p.fgrnhfx,3);
end
[i,j,v]=find(p.fgrnhfx(:,:,nframes));
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
for n=1:nframes
    figure(1)
    subplot(2,2,1)
    mesh(p.fxlong(is:ie,js:je,n),p.fxlat(is:ie,js:je,n),p.fgrnhfx(is:ie,js:je,n)),
    titl('FGRNHFX')
    subplot(2,2,2)
    mesh(p.fxlong(is:ie,js:je,n),p.fxlat(is:ie,js:je,n),p.lfn(is:ie,js:je,n))
    hold on
    contour(p.fxlong(is:ie,js:je,n),p.fxlat(is:ie,js:je,n),p.lfn(is:ie,js:je,n),[0 0])
    % view(0.9,90)
    hold off
    titl('LFN')
    subplot(2,2,3)
    mesh(p.fxlong(is:ie,js:je,n),p.fxlat(is:ie,js:je,n),p.fuel_frac(is:ie,js:je,n))
    titl('FUEL_FRAC')
    % colorbar
    % view(0.9,90)
    subplot(2,2,4)
    d=p.fuel_frac(is:ie,js:je,n);
    if n>1, 
        d=d-p.fuel_frac(is:ie,js:je,n-1);
    end
    mesh(p.fxlong(is:ie,js:je,n),p.fxlat(is:ie,js:je,n),d)
    % colorbar
    titl('FUEL_FRAC diff')
    % view(0.9,90)
    drawnow
    fprintf('frame %i at %s fgrnhfx max %f\n',n,p.times(n,:),p.fgrnhfx_max(n))
    pause(0.3)
end
    function titl(s)  % n and p from outer scope
        title([s,' frame ',num2str(n),' at ',p.times(n,:),' ifire=',num2str(ifire)],...
            'interpreter','none')
    end
end