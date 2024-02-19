function meshmovie(q,name)
% input
% q     structure with fields from wrfout
% name  name of the field to show

v = q.(name);
minv = min(v(:));
maxv = max(v(:));
for i=1:size(q.times,2)
    mesh(q.xlat(:,:,i),q.xlong(:,:,i),v(:,:,i))
    t = char(q.times(:,i)')
    zlim([minv,maxv])
    title(t,interpreter='none')
    drawnow
end