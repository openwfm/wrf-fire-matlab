function saddle_test_n_2d
% dimension
n = [2,2];
fprintf('linear array of %ix%i cells\n',n(1),n(2))

% # fluxes per element
factor = 4;
% # faces = inner faces + the outer face
faces = prod(n)+1; 
% # vertices
vertices = prod(n+1); 
% Euler's formula: v-e+f=2 => e=v+f-2
edges = vertices+faces-2; 
% # constraints for C
c_constraints = edges-n(1)-2*n(2);

% initialize matrices
B=sparse(factor*prod(n));
D=sparse(prod(n),factor*prod(n));
C=sparse(c_constraints,factor*prod(n));
A=sparse(factor*prod(n));
% make them full if small
if prod(n)<200, B=full(B);D=full(D);C=full(C); end
% create block matrices
for i=1:prod(n)
    s=(i-1)*factor+1:i*factor; % span of local dofs
    A(s,s)=diag([1,1,2,2]);
    B(s,s)=[-1,0,0,0;0,1,0,0;0,0,-1,0;0,0,0,1];
    D(i,s)=[1,1,1,1];
end
% define indices of C conditions
inds=index_conditions(n,factor);
% check that they are the same number computed before
if size(inds,1) ~= c_constraints
    error('number of constraints different than indices computed!') 
end
% creating matrix C
for c=1:size(inds,1)
    k=inds(c,inds(c,:)~=0);
    C(c,k) = 1;
end
% initial random wind flux
v0 = rand(factor*prod(n),1);
% solve mass-consistent using saddle problem
saddle_sparse
% show matrix M if it is small
if prod(n)<10,M,end
end

function inds=index_conditions(n,factor)
    % initialize indices matrix
    inds = [];
    % for each row
    for j=1:n(2)
        % x continuity conditions
        t1 = factor*n(1)*(j-1)+2:factor:factor*(n(1)*j-1);
        inds = [inds;[t1',t1'+3]];
        % y continuity conditions (if not last row)
        if j<n(2)
            t2 = factor*n(1)*(j-1)+factor:factor:factor*(n(1)*j);
            inds = [inds;[t2',t2'+factor*(n(1)-1)+3]];
        end
        % ground flux conditions
        t3 = factor*n(1)*j-2;
        inds = [inds;[0,t3]];
    end
end

