% settings
levels = 5;
fp = [2,2,2];

% initialize wind
wp = zeros(levels,3);
for r=1:levels
    % dimension
    n = [5,5,5]*2^(r-1);
    h = [2,2,2]/(2^r);
    factor = 6;
    fprintf('linear array of %ix%ix%i cells\n',n(1),n(2),n(3))
    % creating grid
    X = regular_mesh(n,h,1);
    thx = h(1)*[0:n(1)]'*ones(1,n(2)+1);
    X = add_terrain_to_mesh(X,thx,'shift');
    X = add_terrain_to_mesh(X,'hill','squash',0.1);
    CX = cell(1,3);
    for k = 1:3
        CX{k} = (X{k}(1:end-1,1:end-1,1:end-1)+X{k}(2:end,1:end-1,1:end-1)+X{k}(1:end-1,2:end,1:end-1)+X{k}(1:end-1,1:end-1,2:end)+X{k}(2:end,2:end,1:end-1)+X{k}(2:end,1:end-1,2:end)+X{k}(1:end-1,2:end,2:end)+X{k}(2:end,2:end,2:end))/8;
    end
    xx = CX{1}; yy = CX{2}; zz = CX{3};
    % assembly sparse matrices
    [A,D,E,B,C,v0] = sparse_assembly(X,h);
    % solve mass-consistent using saddle problem
    saddle_sparse
    % transform resulting fluxes into cartesian winds at the centers
    u=E*B*v;
    % create interpolant elements from wind at the center of the cells
    Fx = scatteredInterpolant(CX{1}(:),CX{2}(:),CX{3}(:),u(1:3:end));
    Fy = scatteredInterpolant(CX{1}(:),CX{2}(:),CX{3}(:),u(2:3:end));
    Fz = scatteredInterpolant(CX{1}(:),CX{2}(:),CX{3}(:),u(3:3:end));
    wp(r,1) = Fx(fp); 
    wp(r,2) = Fy(fp); 
    wp(r,3) = Fz(fp);
end
save('test_A','wp');