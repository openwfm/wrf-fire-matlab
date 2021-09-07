function poly_list = make_poly(x,y,rounds)
%x,y lists of points, column vectors
%fills in extra points along a ploygon

p = [x, y];
c = [mean(x) mean(y)];% center point 
v = p-c ; % vectors connecting the central point and the given points 
th = atan2(v(:,2),v(:,1));
[th, idx] = sort(th);   % sort angles 
p = p(idx,:); % sorting the given points
p = [p; p(1,:)];
% figure,plot( p(:,1), p(:,2), '.-r');
% hold on,scatter(p(:,1),p(:,2),'*k'),hold off
% title('First')


%fill in the gaps
for j = 1:rounds
    n = length(p);
    np = [];
    for i = 2:n
        px = 0.5*(p(i,1)+p(i-1,1));
        py = 0.5*(p(i,2)+p(i-1,2));
        np = [np;[px py]];
    end
    p = [p;np];
    %c = [mean(p(:,1)) mean(p(:,2))]; % mean/ central point
    v = p-c ; % vectors connecting the central point and the given points
    th = atan2(v(:,2),v(:,1)); % angle above x axis
    [th, idx] = sort(th);   % sorting the angles
%     figure,plot( p(:,1), p(:,2), '.-r');
%     hold on,scatter(p(:,1),p(:,2),'*k'),hold off
end
%p = unique(p,'rows');
poly_list = p;

    

end % function