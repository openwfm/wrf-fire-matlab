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
    %The slice function requires a meshgrid to work, while X is an ndgrid, so all of the coordinates are transposed.
    %Create a meshgrid with the same array and physical dimensions as the domain arrays
    %CX. Evaluate each grid point at the scattered interpolant
    [CXG_X,CXG_Y,CXG_Z] = meshgrid(linspace(min(X{2}(:)),max(X{2}(:)), n(2)),...
        linspace(min(X{1}(:)),max(X{1}(:)), n(1)), ...
        linspace(min(X{3}(:)),max(X{3}(:)), n(3)));
    
    WX = F1(CXG_Y(:), CXG_X(:), CXG_Z(:));
    WY = F2(CXG_Y(:), CXG_X(:), CXG_Z(:));
    WZ = F3(CXG_Y(:), CXG_X(:), CXG_Z(:));
    
    
    %Convert individual wind vectors (x,y,z) into a full grid with 3-D wind
    %components at each coordinate in the 3-D grid space.
    WX = reshape(WX,n);
    WY = reshape(WY,n);
    WZ = reshape(WZ,n);
    
    
    wind_speed = sqrt(WX.^2 + WY.^2 + WZ.^2);
    
    %Creating contour surfaces and colormaps of the magnitude of the
    %wind-field at 3 different locations
    hsurfaces = slice(CXG_X,CXG_Y,CXG_Z,wind_speed, [ymin, (ymax-ymin)/2, ymax], xmax,zmin);
    set(hsurfaces,'FaceColor','interp','EdgeColor','none')
    if exist('turbo')==2
        colormap turbo
    else
        colormap jet
    end
    hcont = ...
        contourslice(CXG_X,CXG_Y,CXG_Z,wind_speed,[ymin, (ymax-ymin)/2, ymax], xmax,zmin);
    set(hcont,'EdgeColor',[0.7 0.7 0.7],'LineWidth',0.5)
    % % drawing on current figure
    
    colorbar
    axis tight
    title('Streamline Plot with Meshgrid and Wind Color and Contour Map')


    %Note: To use the scattered interpolant function arrays
    %into column-vector format
    if ~exist('scale_t','var') || isempty(scale_t)
        scale_t = 1.1;
    end
    
    if ~exist('t_final','var') || isempty(t_final)
        t_final = scale_t*max(X{1}(:))/mean(W{1}(:));
    end
    
     tspan = [0 t_final];
    
     x0 = 0;
     y0 = X{2}(1,5,1): (n(2) - 10) :X{2}(1,n(2) - 5,1);
     z0 = init_height;


    %plot_mesh_3d(X,[1,nel(1)+1,1,nel(2)+1,2,2])
if ~exist('scale','var')
    scale = 0.9;
end
for i = 1:length(y0)
    for j = 1:length(z0)
        [t,y] = ode23(@(t,y) odefun(t,y, F1, F2, F3), tspan, [x0,y0(i),z0(j)]);
        % Find indices that constrain the streamline to stay inside the
        % domain of the mesh domain.
        count = 1;
        if max(y(:,1)) < max(X{1}(:))
            
            t0 = min(t);
            tf = count*scale_t*max(X{1}(:))/mean(W{1}(:));
            t = [t0, tf];
            [t,y_new] = ode23(@(t,y_new) odefun(t,y_new, F1, F2, F3), tspan, [x0,y0(i),z0(j)]);
            size(y), size(y_new)
            y = cat(1,y,y_new);
            
            count = count+ 1;
        end
        count  = 0;
        ind1 = find(y(:,1)< max(X{1}(:)));
        ind1 = ind1(1: length(ind1) -1);
        ind2 = find(y(:,2)< max(X{2}(:)));
        ind2 = ind2(1: length(ind2) - 1);
        if length(ind1) < length(ind2)
            ind1 = ind1;
        else
            ind1 = ind2;
        end
        
        
        ind1 = find(y(:,1)< max(X{1}(:)));
        ind1 = ind1(1: length(ind1) - 1);
        
        
        hold on
        plot3(y(ind1,1),y(ind1,2), y(ind1,3))
        
        %Plot vectors until second to last index inside the bound to keep
        %vectors in grid
        quiver3(y(ind1,1),y(ind1,2), y(ind1,3), F1(y(ind1,1),y(ind1,2), y(ind1,3)),...
             F2(y(ind1,1),y(ind1,2), y(ind1,3)),F3(y(ind1,1),y(ind1,2), y(ind1,3)), scale)
        
        
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
