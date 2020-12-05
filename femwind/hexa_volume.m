function v=hexa_volume(X)
%   X  columns are corners of hexa in this ordering
%
%   7-----8
%  /     /|
% 5-----6 |
% |       |
% | 3-----4
% |/     /
% 1-----2 

% row = corners of a tetra
tet =  [1 2 3 5
        2 4 3 5
        4 8 6 5
        3 8 7 5
        2 4 6 5
        3 4 8 5];

vis = 0;

v=0;
if vis
    clf,hold off
    for i=1:size(X,2)
        text(X(1,i),X(2,i),X(3,i),num2str(i),'FontSize',16),hold on
    end
    drawnow
end
for i=1:size(tet,1)
    XX = X(:,tet(i,:));
    vv = abs(det([XX;ones(1,4)]))/6;
    if vis
        shp = alphaShape(XX');plot(shp); hold on, drawnow
        i,err_vv = vv - volume(shp),
    end
    v = v+vv;
end % volume of tetra from points ii
if vis, hold off, end