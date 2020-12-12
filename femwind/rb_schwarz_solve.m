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
sub_step=[4,4]
overlap=sub_size-sub_step
% index sets

% divide into intervals in dimensions 1 and 2
sub_nodes=cell(1,2);
for i=1:2
    msub=ceil(n(i)/sub_step(i));
    % sub_nodes{i}=cell(msub,1); 
    for j=1:msub
        from = (j-1)*sub_step(i)+1;
        to = min(from+sub_size(i)-1,n(i));
        sub_nodes{i}{j}=[from:to];
        disp(['direction ',num2str(i),' group ',num2str(j),' indices ',num2str(sub_nodes{i}{j})])
        if to >= n(i)
            nsub(i)=j;
            break
        end
     end
end

check_cover(n(1),sub_nodes{1}(:))
check_cover(n(2),sub_nodes{2}(:))


% global node indices for the subdomains
IX = cell(nsub);
for i1=1:nsub(1)
    for i2=1:nsub(2)
        idx1 = sub_nodes{1}{i1}; % horizontal rectangle
        idx2 = sub_nodes{2}{i2};
        idx3 = 1:n(3); % vertical all
        disp(['subdomain ',num2str([i1,i2]),' nodes ',num2str(idx1),' by ',num2str(idx2)])
        [ix1,ix2,ix3]=ndgrid(idx1,idx2,idx3);
        ix=sub2ind(n,ix1(:),ix2(:),ix3(:));
        IX{i1,i2}=ix;
    end
end

% check indices
check_cover(nn,IX(:))

disp('decomposing subdomain matrices')
R=cell(nsub);  % cholesky factor
P=cell(nsub);  % permutation
for i1=1:nsub(1)
    for i2=1:nsub(2)
        idx1 = sub_nodes{1}{i1}; % horizontal rectangle
        idx2 = sub_nodes{2}{i2};
        idx3 = 1:n(3); % vertical all
        ix = IX{i1,i2};
        [R{i1,i2},ierr,P{i1,i2}]=chol(K(ix,ix),'vector');
        if ierr
            error(['Subdomain matrix is not positive definite, error ',num2str(ierr)])
        end
        p = P{i1,i2}; r = R{i1,i2}; k = K(ix,ix); err(1) = big(k(p,p)-r'*r);
        b = rand(length(p),1); x = r\(r'\b(p)); x(p)=x; err(2) = big(k*x-b);
        fprintf('subdomain %g %g size %g %g %g nodes %g nonzeros %g Cholesky %g error %g %g\n',...
            i1,i2,length(idx1),length(idx2),length(idx3),length(ix),...
            nnz(K(ix,ix)),nnz(R{i1,i2}),err)
    end
end
        
x = zeros(nn,1);
xx = zeros(nn,1);
for it=1:maxit
    for rb1=1:2
        for rb2=1:2
            for i1=rb1:2:nsub(1)
                for i2=rb2:2:nsub(2)
                    % solving horizontal location i1 i2 and vertical line
                    ix = IX{i1,i2};
                    x(ix) = x(ix) + K(ix,ix)\(F(ix) - K(:,ix)'*x);
                    p = P{i1,i2}; % Ksub(p,p)= R'*R
                    res = F(ix) - K(:,ix)'*xx;
                    sol = R{i1,i2}\(R{i1,i2}'\res(p));
                    sol(p) = sol;
                    err = big(res - K(ix,ix)*sol)
                    xx(ix) = xx(ix) + sol;
                    err = big(xx(ix)-x(ix))
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

end

function check_cover(m,x)
% check if cell array x contains indices that cover 1 to m
c=zeros(m,1);
for i=1:length(x)
   c(x{i})=c(x{i})+1;
end
if any(c==0),
    error('subdomains do not cover all nodes')
end
end

