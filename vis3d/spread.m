function spread(p,v)
% spread(p,v)
% display fire spread in variable v of array of structures p
% p(i) needs to have fields v and times
% example:
%   w=wrfout file name
%   p=nc2struct(w,{'FIRE_AREA','Times'},{})
%   spread(p,'fire_area')
    if ~exist('v','var')
        v='fire_area';
    end
    for i=1:length(p),for j=1:size(p(i).fire_area,3)
        contour(p(i).(v)(:,:,j))
        title([v,' ',char(p(i).times(:,j)')],'interpreter','none')
        colorbar
        drawnow
        pause(0.1)
    end
end
