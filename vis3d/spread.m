function spread(p)
% display fire spread in array of structures p
% p(i) needs to have fields fire_area and times
    for i=1:length(p),for j=1:size(p(i).fire_area,3)
        contour(p(i).fire_area(:,:,j))
        title(['FIRE_AREA ',char(p(i).times(:,j)')],'interpreter','none')
        drawnow
    end
end
