function poly_list = make_poly(x,y,rounds)
%x,y lists of points, column vectors

% N = 100 ;
% x = rand(1,N) ;
% y = rand(1,N) ;
p = [x, y]; % coordinates / points 
c = [mean(x) mean(y)]; % mean/ central point 
d = p-c ; % vectors connecting the central point and the given points 
th = atan2(d(:,2),d(:,1)); % angle above x axis
[th, idx] = sort(th);   % sorting the angles 
p = p(idx,:); % sorting the given points
p = [p; p(1,:)]; % add the first at the end to close the polygon 
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
    d = p-c ; % vectors connecting the central point and the given points
    th = atan2(d(:,2),d(:,1)); % angle above x axis
    [th, idx] = sort(th);   % sorting the angles
%     figure,plot( p(:,1), p(:,2), '.-r');
%     hold on,scatter(p(:,1),p(:,2),'*k'),hold off
end
%p = unique(p,'rows');
poly_list = p;

    

end % function