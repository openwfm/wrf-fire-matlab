function ts_at_test
     % TS_AT_TEST tests the TS_AT interpolation function.
    %
    % This function generates a synthetic dataset of variable values defined
    % on a grid with random perturbations of coordinates 
    % and verifies the TS_AT function by checking if the interpolation 
    % at random points preserves linear functions
    %
    % No arguments.
    %
    % No return value.
    %
    % Example:
    %   ts_at_test;  


    % Generate xlon and xlat as uniform with random perturbations
    xlon = repmat(linspace(-180, 180, 100), [100, 1]) + rand(100) - 0.5;
    xlat = repmat(linspace(-90, 90, 100)', [1, 100]) + rand(100) - 0.5;
    
    % Define v as a linear function a*xlon + b*xlat
    a = 2;
    b = 3;
    v=a*xlon + b*xlat;
    v = repmat(v, [1, 1, 10]); % Repeat the same pattern for 10 time steps
    
    % Verify that ts_at returns a*lon + b*lat for several random (lon, lat)
    for i = 1:5
        lon = rand()*360 - 180; % Random longitude between -180 and 180
        lat = rand()*180 - 90;  % Random latitude between -90 and 90
        ts = ts_at(lon, lat, xlon, xlat, v);
        expected_ts = a*lon + b*lat;
        err=abs(ts - expected_ts)
        assert(all(abs(err) < 1e-5), 'Test failed!');
    end
    
    disp('All tests passed.');
end
