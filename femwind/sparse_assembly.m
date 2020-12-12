function [K,F,W]=sparse_assembly(A,X,u0,lambda)
% in: 
%  A  penalty coefficients matrix, size 3x3, s.p.d.
%  X  cell array with 3d arrays x y z coordinates of mesh vertices
%  u0 cell array with 3d arrays initial wind vector at centers in x y z directions
%  lambda scalar field on node
% out:
%  K  stiffness matrix, sparse
%  F  load vector
%  W  gradient of lambda

% do to: now computing all three K,F,W, add a switch to compute only one
% and save time

fprintf('sparse_assembly:')

d = size(X,2);    % dimensions
n = size(X{1});   % mesh size
nn = prod(n);     % total nodes
nz = nn*27;        % estimated nonzeros

% initialize matrices
% KK = sparse([],[],[],nn,nn,nz);
F = zeros(nn,1);
W = {zeros(n-1),zeros(n-1),zeros(n-1)};
nzelem=prod(n-1)*64;     % total nonzeros before assembly
ii=zeros(nzelem,1);
jj=zeros(nzelem,1);
aa=zeros(nzelem,1);
kk=0;

% create matrices
tstart=tic;
Xloc = zeros(3,8);
kglo=zeros(1,8);
m = n-1;
done_last=0;
for i3=1:m(3)
    for i2=1:m(2)
        for i1=1:m(1)  % loop over elements
            % now in element (i1,i2,i3)
            for j3=0:1 % loop over corners of the element
                for j2=0:1
                    for j1=0:1   
                        kloc=sub2ind([2,2,2],j1+1,j2+1,j3+1); % local index
                        % kloc=1+j1+2*(j2+2*j3);   
                        kglo(kloc)=sub2ind(n,i1+j1,i2+j2,i3+j3); % global index
                        for i=1:3
                            Xloc(i,kloc)=X{i}(i1+j1,i2+j2,i3+j3); % node coordinate i
                        end
                    end
                end
            end
            u0loc=[u0{1}(i1,i2,i3),u0{2}(i1,i2,i3),u0{3}(i1,i2,i3)]';
            [Kloc,Floc,Jg]=hexa(A,Xloc,u0loc); % create local matrix and rhs
            F(kglo)=F(kglo)+Floc; % assemble to global rhs
            % KK(kglo,kglo)=KK(kglo,kglo)+Kloc; % assemble to global matrix
            % instead, accumulate contributions to the global matrix
            [ix,jx]=ndgrid(kglo,kglo);
            nzloc=prod(size(Kloc));
            ii(kk+1:kk+nzloc)=ix(:);
            jj(kk+1:kk+nzloc)=jx(:);
            aa(kk+1:kk+nzloc)=Kloc(:);
            kk=kk+nzloc;
            grad = lambda(kglo)'*Jg;  % grad lambda
            grad = grad/A;
            for i=1:3
                W{i}(i1,i2,i3)=grad(i);
            end                         
        end
        % done = 100*((i3-1)*m(2)+i2)/(m(3)*m(2));
        % done = round(done);
        % if done>done_last+5, fprintf(' %g%% ',done), done_last=done; end
    end
end
for i=1:3
    W{i}=u0{i}+W{i};
end
fprintf('\n')
if kk~=nzelem, error('wrong element nonzeros'),end
K = sparse(ii,jj,aa,nn,nn); % err_K=big(K-KK)
nn=prod(n);
nz=nnz(K);
fprintf('stiffness matrix size %g nonzeros %g density %g%%\n',nn,nz,100*nz/nn^2)
% checks
if length(F)>nn, size(F),error('size has grown by an index out of bounds'),end
check_symmetry(K,'K',eps)
K = (K+K')/2; % symmetrize
end
