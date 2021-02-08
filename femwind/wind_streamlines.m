function wind_streamline(X, CX, W, init_height)

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
    %The slice function requires an ndgrid or meshgrid to work. Create a
    %meshgrid with the same array and physical dimensions as the domain arrays
    %CX. Evaluate each grid point at the scattered interpolant
    [CXG_X,CXG_Y,CXG_Z] = meshgrid(linspace(min(X{1}(:)),max(X{1}(:)), n(1)),...
        linspace(min(X{2}(:)),max(X{2}(:)), n(2)), ...
        linspace(min(X{3}(:)),max(X{3}(:)), n(3)));

    WX = F1(CXG_X(:), CXG_Y(:), CXG_Z(:));
    WY = F2(CXG_X(:), CXG_Y(:), CXG_Z(:));
    WZ = F3(CXG_X(:), CXG_Y(:), CXG_Z(:));
    %Convert individual wind vectors (x,y,z) into a full grid with 3-D wind
    %components at each coordinate in the 3-D grid space.
    WX = reshape(WX,n);
    WY = reshape(WY,n);
    WZ = reshape(WZ,n);

    wind_speed = sqrt(WX.^2 + WY.^2 + WZ.^2);
    %Creating countour surfaces and colormaps of the magnitude of the
    %wind-field at 3 different locations
    hsurfaces = slice(CXG_X,CXG_Y,CXG_Z,wind_speed, [xmin, (xmax-xmin)/2, xmax], ymax,zmin);
    set(hsurfaces,'FaceColor','interp','EdgeColor','none')
    colormap turbo
    hcont = ...
        contourslice(CXG_X,CXG_Y,CXG_Z,wind_speed,[xmin,(xmax-xmin)/2,xmax],ymax,zmin);
    set(hcont,'EdgeColor',[0.7 0.7 0.7],'LineWidth',0.5)
    % % drawing on current figure
    
    colorbar
    axis tight
    title('Streamline Plot with Meshgrid and Wind Color and Contour Map')


    %Note: To use the scattered interpolant function arrays
    %into column-vector format

    
    tspan = [0 1000];
    
     x0 = 0;
     y0 = 0: round(max(X{2}(:))/10):max(X{2}(:));
     z0 = init_height;

    %plot_mesh_3d(X,[1,nel(1)+1,1,nel(2)+1,2,2])

for i = 1:length(y0)
    for j = 1:length(z0)
        [t,y] = ode45(@(t,y) odefun(t,y, F1, F2, F3), tspan, [x0,y0(i),z0(j)]);
        % Find indices that constrain the streamline to stay inside the
        % domain of the mesh domain.
        k1 = find(y(:,1)< max(X{1}(:)));
        
        
        hold on
        plot3(y(k1,1),y(k1,2), y(k1,3))
        
        quiver3(y(k1,1),y(k1,2), y(k1,3), F1(y(k1,1),y(k1,2), y(k1,3)),...
             F2(y(k1,1),y(k1,2), y(k1,3)),F3(y(k1,1),y(k1,2), y(k1,3)))
        
        
    end
end

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