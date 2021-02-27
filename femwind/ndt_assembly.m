function K=nd_assembly(A,X,lambda,params);
% in: 
%  A  penalty coefficients matrix, size 3x3, s.p.d.
%  X  cell array with 3d arrays x y z coordinates of mesh vertices
% out:
%  K  stiffness matrix, stored as size (n1,n2,n3,3,3,3)

fprintf('sparse_assembly:')

d = size(X,2);    % dimensions
n = size(X{1});   % mesh size in nodes
nn = prod(n);     % total nodes
kk = 0;

% initialize matrices
F = zeros(nn,1);
K = zeros(n(1),n(2),n(3),3,3,3);

% preallocate matrices
tstart=tic;    % start timer
Xloc = zeros(3,8);
kglo=zeros(1,8);
m = n-1;       % grid size in elements
disp('accumulating global stiffness matrix entries')
for i3=1:m(3)
    for i2=1:m(2)
        for i1=1:m(1)  % loop over elements
            % now in element (i1,i2,i3)
            % build the matrix of node coordinates to pass to hexa
            for j3=0:1 % loop over corners of the element
                for j2=0:1
                    for j1=0:1   
                        jloc=1+j1+2*(j2+2*j3);  % local index of the node in element (i1, i2. i3)
                        for i=1:3
                            Xloc(i,jloc)=X{i}(i1+j1,i2+j2,i3+j3); % node coordinate i
                        end
                    end
                end
            end
            [Kloc,~,~]=hexa(A,Xloc,zeros(3,1)); % compute the local matrix
            % loop over local matrix entries j,k
            for j3=0:1 % 
                for j2=0:1
                    for j1=0:1   
                        jloc=1+j1+2*(j2+2*j3); % local index j
                        for k3=0:1 % loop over corners of the element
                            for k2=0:1
                                for k1=0:1   
                                    kloc=1+k1+2*(k2+2*k3);
                                    % add entry of the local matrix; need
                                    % to offset by 2 because all arrays 
                                    % dimensions in matlab start at 1
                                    K(i1+j1,i2+j2,i3+j3,2+k1-j1,2+k2-j2,2+k3-j3)=...
                                        K(i1+j1,i2+j2,i3+j3,2+k1-j1,2+k2-j2,2+k3-j3)+...
                                        Kloc(jloc,kloc);
                                end
                            end
                        end
                    end
                end
            end                   
        end        
    end
end

end
