function saddle_test_n_2d
% dimension
n = [10,10];
fprintf('linear array of %ix%i cells\n',n(1),n(2))
% grid (regular for the moment)
[x,y] = meshgrid(1:n(1),1:n(2));

%{ 
% too complex
% # faces = inner faces + the outer face
faces = prod(n)+1; 
% # vertices
vertices = prod(n+1); 
% Euler's formula: v-e+f=2 => e=v+f-2
edges = vertices+faces-2;
% # border edges different than the ground
bedges = n(2)+2*n(1);
% # constraints for C
c_constraints = edges-bedges;
%}

% # x continuity constraints
nt1 = (n(1)-1)*n(2);
% # y continuity constraints
nt2 = n(1)*(n(2)-1);
% # ground flux conditions
nt3 = n(2);
c_constraints = nt1+nt2+nt3;

% # fluxes per element
factor = 4;
% initialize matrices
B=sparse(factor*prod(n));
D=sparse(prod(n),factor*prod(n));
C=sparse(c_constraints,factor*prod(n));
A=sparse(factor*prod(n));
E=sparse(size(n,2)*prod(n),factor*prod(n));
v0 = zeros(factor*prod(n),1);
% make them full if small
if prod(n)<200, B=full(B);D=full(D);C=full(C);E=full(E);A=full(A); end
% create block matrices
for i=1:prod(n)
    s=(i-1)*factor+1:i*factor; % span of local dofs
    t=(i-1)*size(n,2)+1:i*size(n,2);
    A(s,s)=diag([1,1,1,1]);
    B(s,s)=[-1,0,0,0;0,1,0,0;0,0,-1,0;0,0,0,1];
    D(i,s)=[1,1,1,1];
    E(t,s)=[-.5,.5,0,0;0,0,-.5,.5];
    v0(s)=[-1,1,0,0]; % initial wind flux
end
% define indices of C conditions
inds=index_conditions_2d(n,factor);
% check that they are the same number computed before
if size(inds,1) ~= c_constraints
    error('number of constraints different than indices computed!') 
end
% creating matrix C
for c=1:size(inds,1)
    k=inds(c,inds(c,:)~=0);
    C(c,k) = 1;
end
u0=E*v0;
figure, quiver(x(:),y(:),u0(1:2:end),u0(2:2:end)), xlabel('x'), ylabel('y'), title('Initial wind')
v0=B*v0; % need initial Cartesian wind instead of flux
% solve mass-consistent using saddle problem
saddle_sparse
% show matrix M if it is small
if prod(n)<10,M,end
% transform fluxes into cartesian winds at the centers
u=E*v;
% plot resulting wind
[x,y] = meshgrid(1:n(1),1:n(2));
figure, quiver(x(:),y(:),u(1:2:end),u(2:2:end)), xlabel('x'), ylabel('y'), title('Mass-consistent solution')
end

function inds=index_conditions_2d(n,factor)
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
    end
    % ground flux conditions
    t3=3:factor:factor*n(1);
    inds = [inds;[zeros(size(t3')),t3']];
end

