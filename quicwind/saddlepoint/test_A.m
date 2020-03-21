function test_A
% settings
levels = 5;
fp = [2,2,2];

% initialize
wp = cell(levels,3);
F = cell(levels,3);
e = zeros(levels,1);
for r=1:levels
    % dimension
    n = [4,4,4]*2^(r-1);
    h = [2,2,2]/(2^r);
    fprintf('linear array of %ix%ix%i cells\n',n(1),n(2),n(3))
    % creating grid
    X = regular_mesh(n,h,1);
    thx = h(1)*[0:n(1)]'*ones(1,n(2)+1);
    X = add_terrain_to_mesh(X,thx,'shift');
   %    X = add_terrain_to_mesh(X,'hill','squash',0.1);
    CX = cell(1,3);
    for k = 1:3
        CX{k} = (X{k}(1:end-1,1:end-1,1:end-1)+X{k}(2:end,1:end-1,1:end-1)+X{k}(1:end-1,2:end,1:end-1)+X{k}(1:end-1,1:end-1,2:end)+X{k}(2:end,2:end,1:end-1)+X{k}(2:end,1:end-1,2:end)+X{k}(1:end-1,2:end,2:end)+X{k}(2:end,2:end,2:end))/8;
    end
    % assembly sparse matrices
    [A,D,E,B,C,v0] = sparse_assembly(X,h);
    % solve mass-consistent using saddle problem
    saddle_sparse
    % transform resulting fluxes into cartesian winds at the centers
    u=E*B*v;
    % save errors
    e(r) = err;
    % create interpolant elements from wind at the center of the cells
    F{r,1} = scatteredInterpolant(CX{1}(:),CX{2}(:),CX{3}(:),u(1:3:end));
    F{r,2} = scatteredInterpolant(CX{1}(:),CX{2}(:),CX{3}(:),u(2:3:end));
    F{r,3} = scatteredInterpolant(CX{1}(:),CX{2}(:),CX{3}(:),u(3:3:end));
    wp{r,1} = F{r,1}(fp); 
    wp{r,2} = F{r,2}(fp); 
    wp{r,3} = F{r,3}(fp);
end
d.F = F; d.wp = wp; d.err = e;
save('test_A_shift_only_5','-struct','d');
end