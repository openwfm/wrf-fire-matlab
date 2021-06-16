% User options
wrfout='./wrfout';
s=2;               % sample in all directions 
wlevel=1;          % wind level to plot field
xl=[1600,3600];    % x limits in meters
yl=[600,2600];     % y limits in meters

% Plotting limits chosen
figure, 
h=surf(X{1}(:,:,1),X{2}(:,:,1),X{3}(:,:,1)); 
set(h,'EdgeColor','None'); 
hold on; 
plot3([xl(1),xl(1),xl(2),xl(2),xl(1)],[yl(1),yl(2),yl(2),yl(1),yl(1)],[10000,10000,10000,10000,10000],'k'); 
view(2); 
xlabel('X (m)'); 
ylabel('Y (m)');

% Calculating mesh-center coordinates
[CX,CH] = center_mesh(X); 
fix=find(abs(X{1}(:,1,1)-xl(1)) == min(abs(X{1}(:,1,1)-xl(1))));
fiy=find(abs(X{2}(1,:,1)-yl(1)) == min(abs(X{2}(1,:,1)-yl(1))));
ffx=find(abs(X{1}(:,1,1)-xl(2)) == min(abs(X{1}(:,1,1)-xl(2))));
ffy=find(abs(X{2}(1,:,1)-yl(2)) == min(abs(X{2}(1,:,1)-yl(2))));

% Plotting initial wind and massconsistent wind slice on terrain
az=-130;          % azimuth camera angle
el=15;            % elevation camera angle
sr=ceil((ffx+fix)/2);       % x-slice
% initial wind field
xx=squeeze(CX{1}(sr,fiy:s:ffy,1:s:end)); % x-array for quiver
yy=squeeze(CX{2}(sr,fiy:s:ffy,1:s:end)); % y-array for quiver
zz=squeeze(CX{3}(sr,fiy:s:ffy,1:s:end)); % z-array for quiver
uu=squeeze(u0{1}(sr,fiy:s:ffy,1:s:end)); % u-array for quiver
vv=squeeze(u0{2}(sr,fiy:s:ffy,1:s:end)); % v-array for quiver
ww=squeeze(u0{3}(sr,fiy:s:ffy,1:s:end)); % w-array for quiver
figure,
h=mesh(squeeze(X{1}(:,:,1)),squeeze(X{2}(:,:,1)),squeeze(X{3}(:,:,1))); % plotting terrain 
set(h,'FaceAlpha',.4); hold on; % set transparency and continue plotting
quiver3(xx,yy,zz,uu,vv,ww,'b'); % plotting wind 
view(az,el); xlabel('X (m)'); ylabel('Y (m)'); zlabel('Elevation (m)'); xlim(xl); ylim(yl); % setting view, labels, and axes limits
% mass-consistent wind field
uu=squeeze(W{1}(sr,fiy:s:ffy,1:s:end)); % u-array for quiver
vv=squeeze(W{2}(sr,fiy:s:ffy,1:s:end)); % v-array for quiver
ww=squeeze(W{3}(sr,fiy:s:ffy,1:s:end)); % w-array for quiver
figure, 
h=mesh(squeeze(X{1}(:,:,1)),squeeze(X{2}(:,:,1)),squeeze(X{3}(:,:,1))); % plotting terrain  
set(h,'FaceAlpha',.4); hold on; % set transparency and continue plotting
quiver3(xx,yy,zz,uu,vv,ww,'b'); % plotting wind 
view(az,el); xlabel('X (m)'); ylabel('Y (m)'); zlabel('Elevation (m)'); xlim(xl); ylim(yl); % setting view, labels, and axes limits

% Plotting surface initial wind and massconsistent wind on terrain
az=25;             % azimut angle of camera
el=54;             % elevation angle of camera
% initial wind field
xx=squeeze(CX{1}(fix:s:ffx,fiy:s:ffy,wlevel)); % x-array for quiver
yy=squeeze(CX{2}(fix:s:ffx,fiy:s:ffy,wlevel)); % y-array for quiver
zz=squeeze(CX{3}(fix:s:ffx,fiy:s:ffy,wlevel)); % z-array for quiver
uu=squeeze(u0{1}(fix:s:ffx,fiy:s:ffy,wlevel)); % u-array for quiver
vv=squeeze(u0{2}(fix:s:ffx,fiy:s:ffy,wlevel)); % v-array for quiver
ww=squeeze(u0{3}(fix:s:ffx,fiy:s:ffy,wlevel)); % w-array for quiver
figure,
h=mesh(squeeze(X{1}(:,:,1)),squeeze(X{2}(:,:,1)),squeeze(X{3}(:,:,1))); % plotting terrain 
set(h,'FaceAlpha',.4); hold on; % set transparency and continue plotting
quiver3(xx,yy,zz,uu,vv,ww,'b'); % plotting wind 
view(az,el); xlabel('X (m)'); ylabel('Y (m)'); zlabel('Elevation (m)'); xlim(xl); ylim(yl); % setting view, labels, and axes limits
% mass-consistent wind field
xx=squeeze(CX{1}(fix:s:ffx,fiy:s:ffy,wlevel)); % x-array for quiver
yy=squeeze(CX{2}(fix:s:ffx,fiy:s:ffy,wlevel)); % y-array for quiver
zz=squeeze(CX{3}(fix:s:ffx,fiy:s:ffy,wlevel)); % z-array for quiver
uu=squeeze(W{1}(fix:s:ffx,fiy:s:ffy,wlevel)); % u-array for quiver
vv=squeeze(W{2}(fix:s:ffx,fiy:s:ffy,wlevel)); % v-array for quiver
ww=squeeze(W{3}(fix:s:ffx,fiy:s:ffy,wlevel)); % w-array for quiver
figure,
h=mesh(squeeze(X{1}(:,:,1)),squeeze(X{2}(:,:,1)),squeeze(X{3}(:,:,1))); % plotting terrain 
set(h,'FaceAlpha',.4); hold on; % set transparency and continue plotting
quiver3(xx,yy,zz,uu,vv,ww,'b'); % plotting wind 
view(az,el); xlabel('X (m)'); ylabel('Y (m)'); zlabel('Elevation (m)'); xlim(xl); ylim(yl); % setting view, labels, and axes limits
