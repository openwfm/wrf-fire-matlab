function ts = ts_at(lon, lat, xlon, xlat, v)
    % TS_AT interpolates time series of a variable from a grid to a point.
    %
    % This function performs interpolation of a variable defined on a 2D grid
    % to a specific geographic point (longitude, latitude) over a series of 
    % time steps.
    %
    % Arguments:
    %   lon (1D array) - The longitude of the location to interpolate to.
    %   lat (double) - The latitude of the location to interpolate to.
    %   xlon (2D double array) - The longitudes of the grid points.
    %   xlat (2D double array) - The latitudes of the grid points.
    %   v (3D double array) - The values on the grid. The third dimension 
    %                         is time.
    %
    % Returns:
    %   ts (1D double array) - The interpolated time series of the variable 
    %                           at the specified location.
    %
    % Example:
    %   % Given a set of grid points (xlon, xlat) and variable values v(time)
    %   lon = -120.5; lat = 34.5;
    %   ts_values = ts_at(lon, lat, grid_lon, grid_lat, variable_values);

    % Check input dimensions
    if ndims(v) ~= 3
        error('The variable v must be a 3D array where the last dimension is time.');
    end
    
    % Find the four nearest grid points
    [numRows, numCols] = size(xlon);
    [~, idx] = min(abs(xlon(:) - lon) + abs(xlat(:) - lat));
    [iRow, iCol] = ind2sub([numRows, numCols], idx);
    
    % Define the square subgrid indices
    iRow = max(iRow-1, 1):min(iRow+1, numRows);
    iCol = max(iCol-1, 1):min(iCol+1, numCols);
    
    % Extract the subgrid coordinates
    subXlon = xlon(iRow, iCol);
    subXlat = xlat(iRow, iCol);
    
    % Initialize the output time series
    numTimeSteps = size(v, 3);
    ts = zeros(numTimeSteps, 1);
    
    % Loop over each time step
    for t = 1:numTimeSteps
        % Extract the values for the current time step
        subV = v(iRow, iCol, t);
        
        % Flatten the subgrid values to column vectors
        subV = subV(:);
        
        % Create an interpolant for the current time step
        F = scatteredInterpolant(subXlon(:), subXlat(:), subV);
        
        % Perform the interpolation for the current time step
        ts(t) = F(lon, lat);
    end
end
