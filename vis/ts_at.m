function ts = ts_at(lon, lat, xlon, xlat, v)
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
    
    % Create interpolant with coordinates
    F = scatteredInterpolant(subXlon(:), subXlat(:), []);
    
    % Loop over each time step
    for t = 1:numTimeSteps
        % Extract the values for the current time step
        subV = v(iRow, iCol, t);
        
        % Assign the values to the interpolant and interpolate
        F.Values = subV(:);
        ts(t) = F(lon, lat);
    end
end
