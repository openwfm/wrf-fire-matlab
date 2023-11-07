function ts_at_test
    % Test the ts_at function
    
    % Generate xlon and xlat as uniform with random perturbations
    xlon = repmat(linspace(-180, 180, 100), [100, 1]) + rand(100) - 0.5;
    xlat = repmat(linspace(-90, 90, 100)', [1, 100]) + rand(100, 1) - 0.5;
    
    % Define v as a linear function a*xlon + b*xlat
    a = 2;
    b = 3;
    v = bsxfun(@plus, bsxfun(@times, a, xlon), bsxfun(@times, b, xlat));
    v = repmat(v, [1, 1, 10]); % Repeat the same pattern for 10 time steps
    
    % Verify that ts_at returns a*lon + b*lat for several random (lon, lat)
    for i = 1:5
        lon = rand()*360 - 180; % Random longitude between -180 and 180
        lat = rand()*180 - 90;  % Random latitude between -90 and 90
        ts = ts_at(lon, lat, xlon, xlat, v);
        expected_ts = a*lon + b*lat;
        assert(all(abs(ts - expected_ts) < 1e-5), 'Test failed!');
    end
    
    disp('All tests passed.');
end
