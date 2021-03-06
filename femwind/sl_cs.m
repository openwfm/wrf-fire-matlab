function stream_cs(X, CX, W, init_z, n_lines, end_z, y0)
%X - Wind mesh grid
%CX Center coordinates of each grid element
%W Wind Array Used For Streamlines
%Parameters set by user
% init_z: Initial height for plotting first streamline
% n_lines: number of streamlines
% end_z: 
%W Final wind, or gradient matrix at center points of the grid elements
%disp('wind_streamlines not done yet')
%Extracting min and max values from mesh grid to create a valid n-d grid 
%for streamlines function
% Most work goes into creating a grid space with
    xmin = min(X{1}(:));
    xmax = max(X{1}(:));
    ymin = min(X{2}(:));
    ymax = max(X{2}(:));
    zmin = min(X{3}(:));
    zmax = max(X{3}(:));
    %Creating nd-meshgrid
    [n(1),n(2),n(3)]=size(X{1});

    level = 1:n(3);
    [F1, F2, F3] = wind_interp(W,CX);
    %The slice function requires a meshgrid to work, while X is an ndgrid, so all of the coordinates are transposed.
    %Create a meshgrid with the same array and physical dimensions as the domain arrays
    %CX. Evaluate each grid point at the scattered interpolant
%     [CXG_X,CXG_Y,CXG_Z] = meshgrid(linspace(min(X{2}(:)),max(X{2}(:)), n(2)),...
%         linspace(min(X{1}(:)),max(X{1}(:)), n(1)), ...
%         linspace(min(X{3}(:)),max(X{3}(:)), n(3)));
%     
%     WX = F1(CXG_Y(:), CXG_X(:), CXG_Z(:));
%     WY = F2(CXG_Y(:), CXG_X(:), CXG_Z(:));
%     WZ = F3(CXG_Y(:), CXG_X(:), CXG_Z(:));
%     
%     
%     %Convert individual wind vectors (x,y,z) into a full grid with 3-D wind
%     %components at each coordinate in the 3-D grid space.
%     WX = reshape(WX,n);
%     WY = reshape(WY,n);
%     WZ = reshape(WZ,n);
%     
%     
%     wind_speed = sqrt(WX.^2 + WY.^2 + WZ.^2);
%     
%     %Creating contour surfaces and colormaps of the magnitude of the
%     %wind-field at 3 different locations
%     hsurfaces = slice(CXG_X,CXG_Y,CXG_Z,wind_speed, [ymin, (ymax-ymin)/2, ymax], xmax,zmin);
%     set(hsurfaces,'FaceColor','interp','EdgeColor','none')
%     colormap turbo
%     hcont = ...
%         contourslice(CXG_X,CXG_Y,CXG_Z,wind_speed,[ymin, (ymax-ymin)/2, ymax], xmax,zmin);
%     set(hcont,'EdgeColor',[0.7 0.7 0.7],'LineWidth',0.5)
%     % % drawing on current figure
%     
%     colorbar
%     axis tight
%     title('Streamline Plot with Meshgrid and Wind Color and Contour Map')


    %Note: To use the scattered interpolant function arrays
    %into column-vector format
     
     tspan = [0 100000];
     x0 = 0;
     y0 = y0; %(max(X{2}(:)) + min(X{2}(:)))/2 + 1;
     %z0 = init_z:n_lines: end_z;
     z0 = 2:25: max(X{3}(:));


    %plot_mesh_3d(X,[1,nel(1)+1,1,nel(2)+1,2,2])
if ~exist('scale','var')
    scale = 0.9;
end
for i = 1:length(y0)
    for j = 1:length(z0)
        [t,y] = ode45(@(t,y) odefun(t,y, F1, F2, F3), tspan, [x0,y0(i),z0(j)]);
        % Find indices that constrain the streamline to stay inside the
        % domain of the mesh domain.
        k1 = find(y(:,1)< max(X{1}(:)));
        k1 = k1(1: length(k1));
        
        hold on
        plot(y(k1,1), y(k1,3))
        
        %Plot vectors until second to last index inside the bound to keep
        %vectors in grid
%         quiver(y(k1,1), y(k1,3),...
%              F1(y(k1,1),y(k1,2), y(k1,3)),F3(y(k1,1),y(k1,2), y(k1,3)))
        
    title(['XZ-Planar Cross Sectional Profile of Streamlines Positioned at Y=', num2str(y0)])
    xlabel('X')
    ylabel('Z')
    end
end
axis tight


%%
    function [F1,F2,F3] = wind_interp(W, CX)
        
        F1 = scatteredInterpolant(CX{1}(:), CX{2}(:), CX{3}(:),W{1}(:));
        F2 = scatteredInterpolant(CX{1}(:), CX{2}(:), CX{3}(:),W{2}(:));
        F3 = scatteredInterpolant(CX{1}(:), CX{2}(:), CX{3}(:),W{3}(:));
    end



    function dXdt = odefun(t,y, F1, F2, F3)
        dXdt = zeros(3,1);
        dXdt(1) = F1(y(1),y(2),y(3));
        dXdt(2) = F2(y(1),y(2),y(3));
        dXdt(3) = F3(y(1),y(2),y(3));
    end
end 