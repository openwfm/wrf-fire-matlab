function [XH,WH]=wind_at_h(X,CX,W,m,bbox)
% [XH,WH]=wind_at_h(X,CX,W,m,bbox)
% interpolate and plot wind at model independent mesh size m(1) m(2) m(3)
% in a bbox above terrain

% create query mesh
[XH{1},XH{2},XH{3}]=ndgrid(linspace(bbox(1),bbox(2),m(1)),...
                           linspace(bbox(3),bbox(4),m(2)),...
                           linspace(bbox(5),bbox(6),m(3)));
% interpolate terrain from nodes
terrain = griddata(X{1}(:,:,1),X{2}(:,:,1),X{3}(:,:,1),XH{1},XH{2});
for i=1:m(3),
    XH{3}(:,:,i)=XH{3}(:,:,i)+terrain;
end
% interpolate wind components from cell centers
for i=1:3
    WH{i}=griddata(CX{1},CX{2},CX{3},W{i},XH{1},XH{2},XH{3});
end
plot_wind_3d(XH,WH)
end
