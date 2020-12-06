function C = sparse_C(n)
% C = sparse_C(n)
% Create the wind face continuity constraint matrix for a mesh of a given size
% in:
%     n   vector size 3, mesh sizes in x y z directions
% out
%     C   sparse matrix, one row for each cell face including on the ground

% dimension and factor
d = size(n,2);
factor = 2*d;
% # ground flux conditions
ncg = n(1)*n(2);
% # continuity constraints
ncx = (n(1)-1)*n(2)*n(3); 
ncy = n(1)*(n(2)-1)*n(3); 
ncz = n(1)*n(2)*(n(3)-1);
% total number of constraints
c_constraints = ncg + ncx + ncy + ncz;
% initialize arrays
cg = sparse(ncg,factor*prod(n)); cx = sparse(ncx,factor*prod(n));
cy = sparse(ncy,factor*prod(n)); cz = sparse(ncz,factor*prod(n));
for e=1:prod(n)
    [xi,yi,zi]=ind2sub(n,e); % 3D indexing from 1D
    s=(e-1)*factor+1:e*factor; % span of local dofs
    % continuity and ground flux conditions
    [cx,cy,cz,cg]=c_conditions_3d(n,xi,yi,zi,s,cx,cy,cz,cg);
end
% full continuity operator
C = [cx;cy;cz;cg];
% check number of continuity constraints
if size(C,1) ~= c_constraints
    error('number of constraints different than indices computed!') 
end
end