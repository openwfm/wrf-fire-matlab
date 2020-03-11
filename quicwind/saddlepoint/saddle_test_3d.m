% number of cells and spacing
n = [5,5,5];
h = [1,1,1];
fprintf('linear array of %ix%ix%i cells\n',n(1),n(2),n(3))
% dimensions
d = size(n,2);
% number of dof per cell
factor = 2*d;

% creating grid
X = regular_mesh(n,h,1);
X = add_terrain_to_mesh(X,'hill','squash',0.1);
CX = center_mesh(X);
xx = CX{1}; yy = CX{2}; zz = CX{3};

% initial wind
v0 = zeros(factor*prod(n),1);
for e=1:prod(n)
    s=(e-1)*factor+1:e*factor; % span of local dofs
    % initial wind
    v0(s)=[1,1,0,0,0,0]';
end

% create continuity and ground flux conditions
C = sparse_C(n);

% static cell matrices
E = [.5,.5,0,0,0,0;
    0,0,.5,.5,0,0;
    0,0,0,0,.5,.5]; % matrix of interpolated winds at the center of the cells
D = ones(1,factor); % diverge operator

% moduli for penalization
moduli = [1e3,1e3,1];
% initialize resulting wind
Pf = sparse(factor*prod(n),factor*prod(n));
v = zeros(factor*prod(n),1);
% solve mass-consistent using 3D indexing
for k=1:n(3)
    for j=1:n(2)
        for i=1:n(1)
            ind = sub2ind(n,i,j,k);
            s=(ind-1)*factor+1:ind*factor;
            [A,B] = cell_A_B(h(1),h(2),X{3},i,j,k,moduli);
            Pf(s,s) = cell_P(A,B,D);
        end
    end
end