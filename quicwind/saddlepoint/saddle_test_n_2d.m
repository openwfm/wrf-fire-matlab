function saddle_test_n_2d
n = [2,2];

ffactor = 1;
efactor = 4;
faces = ffactor*(prod(n)+1); % inner faces + the outer face
vertices = prod(n+1); 
% Euler's formula: v-e+f=2 => e=v+f-2
edges = vertices+faces-2; 
c_constraints = edges-n(1)-2*n(2);

fprintf('linear array of %ix%i cells\n',n(1),n(2))
B=sparse(efactor*prod(n));
D=sparse(prod(n),efactor*prod(n));
C=sparse(c_constraints,efactor*prod(n));
A=sparse(efactor*prod(n));
if prod(n)<200, B=full(B);D=full(D);C=full(C); end
for i=1:prod(n)
    s=(i-1)*efactor+1:i*efactor; % span of local dofs
    A(s,s)=diag([1,1,2,2]);
    B(s,s)=[-1,0,0,0;0,1,0,0;0,0,-1,0;0,0,0,1];
    D(i,s)=[1,1,1,1];
end
inds=index_conditions(n,efactor);
if size(inds,1) ~= c_constraints
    error('number of constraints different than indices computed!') 
end
for c=1:size(inds,1)
    k=inds(c,inds(c,:)~=0);
    C(c,k) = 1;
end
v0 = rand(efactor*prod(n),1);
saddle_sparse
if prod(n)<10,M,end

function inds=index_conditions(n,efactor)
    inds = [];
    for j=1:n(2)
        t1 = efactor*n(1)*(j-1)+2:efactor:efactor*(n(1)*j-1);
        inds = [inds;[t1',t1'+3]];
        if j<n(2)
            t2 = efactor*n(1)*(j-1)+efactor:efactor:efactor*(n(1)*j);
            inds = [inds;[t2',t2'+efactor*(n(1)-1)+3]];
        end
        t3 = efactor*n(1)*j-2;
        inds = [inds;[0,t3]];
    end
end
end

