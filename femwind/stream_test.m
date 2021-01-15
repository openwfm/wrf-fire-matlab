% load wind
disp('wind_streamlines not done yet')
% x = X{1}(1:132,1:132,1:16);
% y = X{2}(1:132,1:132,1:16);
% z = X{3}(1:132,1:132,1:16);
%POssibl


%Extracting min and max values from mesh grid to create a valid n-d grid 
%for streamlines function
xmin = min(CX{1}(:)); 
xmax = max(CX{1}(:)); 
ymin = min(CX{2}(:));
ymax = max(CX{2}(:));
zmin = min(CX{3}(:));
zmax = max(CX{3}(:));
%Creating nd-meshgrid
[n(1),n(2),n(3)]=size(CX{1});
%[x,y,z] = meshgrid(xmin:5:xmax, ymin:5:ymax, 1:1:n(3));
level = 1:n(3);

wind_speed = sqrt(W{1}.^2 + W{2}.^2 + W{3}.^2);
%wind_speed = sqrt(u.^2 + v.^2 + w.^2);
%Issue with x,y,z not being in an NDGRID form?
hsurfaces = slice(CX{1}(:,:, level),CX{2}(:,:, level),CX{3}(:,:, level), wind_speed,ctour,ymax,zmin);
set(hsurfaces,'FaceColor','interp','EdgeColor','none')
colormap turbo
hcont = ...
contourslice(x,y,z,wind_speed,[xmin,400,xmax],ymax,zmin);
set(hcont,'EdgeColor',[0.7 0.7 0.7],'LineWidth',0.5)
% % drawing on current figure
% title('wind_streamlines not done yet')
 [sx,sy,sz] = meshgrid(0,0:25:650,1);
 hlines = streamline(CX{1}(:,:,level),CX{2}(:,:,level),CX{3}(:,:,level),...
   W{1}(:,:,level), W{2}(:,:,level), W{3}(:,:,level),sx,sy,sz);
 set(hlines,'LineWidth',2,'Color','r')
% 
view(3)
daspect([40,40,1])
axis tight
% title('wind_streamlines not done yet')  


%%
%Note: To use the scattered interpolant function arrays 
%into column-vector format

[F1, F2, F3] = wind_interp(W,CX);
tspan = [0 600];
% x0 = (max(XC{1}(:)) - min(XC{1}(:)))/2;
% y0 = (max(XC{2}(:)) - min(XC{2}(:)))/2;
% z0 = (max(XC{3}(:)) - min(XC{3}(:)))/2; %initial values for coordinates not being evolved in time

x0 = 1;
y0 = 0:5:650;
z0 = 5;
init_pos = [x0,y0,z0];


plot_mesh_3d(X,[1,nel(1)+1,1,nel(2)+1,1,1])

for i = 1:length(y0)
[t,y] = ode45(@(t,y) odefun(t,y, F1, F2, F3), tspan, [x0,y0(i),z0])
% [t,x] = ode45(@(t,x) odesolv1(t,x, y0, z0, F1), tspan, x0)
% [t2,y] = ode45(@(t,y) odesolv2(t,x0, y, z0, F2), tspan, y0)
% [t3,z] = ode45(@(t,z) odesolv3(t,x0, y0, z, F3), tspan, z0)
% figure(8)
% height=10;
% bbox=[min(X{1}(:)),max(X{1}(:)),...
% min(X{2}(:)),max(X{2}(:)),...
% height,height];
% [XH,WH]=wind_at_h(X,CX,W,[20,20,1],bbox);
% plot_wind_3d(XH,WH)
% hold on

hold on
plot3(y(:,1),y(:,2), y(:,3))
end
%plot3(x,y,z)
function [F1,F2,F3] = wind_interp(W, CX)

    F1 = scatteredInterpolant(CX{1}(:), CX{2}(:), CX{3}(:),W{1}(:));
    F2 = scatteredInterpolant(CX{1}(:), CX{2}(:), CX{3}(:),W{2}(:));
    F3 = scatteredInterpolant(CX{1}(:), CX{2}(:), CX{3}(:),W{3}(:));

%      F1 = scatteredInterpolant(CX{1}(:,:,1:2), CX{2}(:,:,1:2), CX{3}(:,:,1:2),W{1}(:,:,1:2));
%      F2 = scatteredInterpolant(CX{1}(:,:,1:2), CX{2}(:,:,1:2), CX{3}(:,:,1:2),W{2}(:,:,1:2));
%      F3 = scatteredInterpolant(CX{1}(:,:,1:2), CX{2}(:,:,1:2), CX{3}(:,:,1:2),W{3}(:,:,1:2));
end


function dXdt = odefun(t,y, F1, F2, F3)
    dXdt = zeros(3,1);
    dXdt(1) = F1(y(1),y(2),y(3));
    dXdt(2) = F2(y(1),y(2),y(3));
    dXdt(3) = F3(y(1),y(2),y(3));
end
