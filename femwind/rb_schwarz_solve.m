function x=rb_schwarz_solve(K,F,X)
% x=rb_line_gs_solve(K,F,X)
disp('solver red-black horizonal schwarz')
n = size(X{1})
nn = size(F,1)
if nn~=prod(n)
    error('rb_line_gs_solve: inconsistent sizes')
end

% parameters
maxit=1000
tol = 1e-5*big(F)/big(K)
sub_size=[6,6]
sub_step=[5,5]
overlap=sub_size-sub_step
% index sets

% divide indices
sub_nodes=cell(1,2);
for i=1:2
    nsub(i)=ceil(n(i)/sub_step(i));
    sub_nodes{i}=cell(nsub(i),1); 
    for j=1:nsub(i)
        ix = (j-1)*sub_step(i);
        from = ix+1;
        to = min(ix+sub_size(i),n(i));
        sub_nodes{i}{j}=[from:to];
        disp(['direction ',num2str(i),' group ',num2str(j),' indices ',num2str(sub_nodes{i}{j})])
    end
end

% decompose subdomain matrices
R=cell(nsub);  % cholesky factor
P=cell(nsub);  % permutation
ix=cell(nsub); % indices of selection
for i1=1:nsub(1)
    for i2=1:nsub(2)
        idx1 = sub_nodes{1}{i1}; % horizontal rectangle
        idx2 = sub_nodes{2}{i2};
        idx3 = 1:n(3); % vertical all
        disp(['subdomain ',num2str([i1,i2]),' nodes ',num2str(idx1),' by ',num2str(idx2)])
        [ix1,ix2,ix3]=ndgrid(idx1,idx2,idx3);
        ix=sub2ind(n,ix1(:),ix2(:),ix3(:));
        [R{i1,i2},FLAG,P{i1,i2}]=chol(K(ix,ix));
        if FLAG,
            error('Subdomain matrix is not positive definite')
        end
        fprintf('subdomain %g %g size %g %g %g nodes %g nonzeros %g Cholesky %g\n',...
            i1,i2,length(idx1),length(idx2),length(idx3),length(ix),...
            nnz(K(ix,ix)),nnz(R{i1,i2}))
    end
end
        
        

x = zeros(nn,1);
colx=1:n(3);
onex=ones(1,n(3));
for it=1:maxit
    for rb1=1:2
        for rb2=1:2
            for i1=rb1:2:n(1)
                for i2=rb2:2:n(2)
                    % solving horizontal location i1 i2 and vertical line
                    ix = sub2ind(n,i1*onex,i2*onex,colx); 
                    x(ix) = x(ix) - K(ix,ix)\(K(:,ix)'*x - F(ix));
                end
            end
        end
    end
    res= norm(K*x-F);
    fprintf('iteration %g residual %g tolerance %g\n',it,res,tol)
    if res<tol,
        break
    end
end