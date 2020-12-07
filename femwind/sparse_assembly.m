function [K,F]=sparse_assembly_fem(A,X,u0)
% in: 
%  A  penalty coefficients matrix, size 3x3, s.p.d.
%  X  cell array with 3d arrays x y z coordinates of mesh vertices
%  u0 cell array with 3d arrays initial wind vector at centers in x y z directions

d = size(X,2);    % dimensions
n = size(X{1});   % mesh size
nn = prod(n);     % total nodes
nz = nn*8;        % estimated nonzeros

% initialize matrices
K = sparse([],[],[],nn,nn,nz);
F = zeros(nn,1);

% create matrices
tstart=tic;
Xloc = zeros(3,8);
kglo=zeros(1,8);
for i1=1:n(1)-1  % loop over elements
    for i2=1:n(2)-1
        for i3=1:n(3)-1
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
            [Kloc,Floc]=hexa(A,Xloc,u0loc); % create local matrix and rhs
            K(kglo,kglo)=K(kglo,kglo)+Kloc; % assemble to global matrix
            F(kglo)=F(kglo)+Floc; % assemble to global rhs
        end
    end
end
% checks
if length(F)>nn, size(F),error('size has grown by an index out of bounds'),end
check_symmetry(K,'K',eps)
end
