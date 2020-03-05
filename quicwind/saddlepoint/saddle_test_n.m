n=5
fprintf('linear array of %i cells\n',n)
B=sparse(2*n);
D=sparse(n,2*n);
C=sparse(n,2*n);
A=speye(2*n);
if n<200, B=full(B);D=full(D);C=full(C); end
for i=1:n
    s=(i-1)*2+1:i*2;    % span of local dofs
    A(s,s)=diag([1,2]);
    B(s,s)=[-1,0;0,1];
    D(i,s)=[1,1];
    C(i,s(end))=1;
    if i<n,
        C(i,s(end)+1)=1; 
    else
        disp('flux boundary condition on last cell')
    end
end
v0 = rand(2*n,1);
saddle_sparse
if n<10,M,end
