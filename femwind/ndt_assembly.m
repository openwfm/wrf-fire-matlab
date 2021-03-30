function K=ndt_assembly(A,X,u0,lambda,params,m);
% in: 
%  A  penalty coefficients matrix, size 3x3, s.p.d.
%  X  cell array with 3d arrays x y z coordinates of mesh vertices
%  u0, lambda, params  not used
%  storage scheme, 14 or 27
% out:
%  K  stiffness matrix, stored as size (n1,n2,n3,3,3,3)

d = size(X,2);    % dimensions
n = size(X{1});   % mesh size in nodes
nn = prod(n);     % total nodes

t = ndt_storage_table(m); 

% initialize matrices
F = zeros(nn,1);
K = zeros(n(1),n(2),n(3),m);
Xloc = zeros(3,8); 
for ie3=1:n(3)-1  % loop over elements
    for ie2=1:n(2)-1
        for ie1=1:n(1)-1  
            % now in element (ie1,ie2,ie3)
            % build the matrix of node coordinates to pass to hexa
            for ic3=0:1 % loop over corners of the element
                for ic2=0:1
                    for ic1=0:1   
                        iloc=1+ic1+2*(ic2+2*ic3);  % local index of the node in the element
                        for i=1:3
                            Xloc(i,iloc)=X{i}(ie1+ic1,ie2+ic2,ie3+ic3); % node coordinate i
                        end
                    end
                end
            end
            [Kloc,~,~]=hexa(A,Xloc,zeros(3,1)); % compute the local stiffness matrix
            % loop over element corners ic, kc
            for ic3=0:1 % 
                for ic2=0:1
                    for ic1=0:1   
                        for kc3=0:1 
                            for kc2=0:1
                                for kc1=0:1   
                                    iloc=1+ic1+2*(ic2+2*ic3); % index in the local element matrix
                                    kloc=1+kc1+2*(kc2+2*kc3); % index in the local element matrix
                                    % global index triple of node i
                                    i1=ie1+ic1;
                                    i2=ie2+ic2;
                                    i3=ie3+ic3;
                                    % relative position of k vs. i
                                    j1 = kc1-ic1;
                                    j2 = kc2-ic2;
                                    j3 = kc3-ic3;
                                    % index triple of row m of K where the entry (i,k) is stored
                                    % in fortran we won't have the 2+ because
                                    % the array t will be indexed -1:1
                                    % storing in this row only 
                                    m3 = i3+t(3,2+j1,2+j2,2+j3);
                                    m2 = i2+t(2,2+j1,2+j2,2+j3);
                                    m1 = i1+t(1,2+j1,2+j2,2+j3);
                                    % index of the matrix entry (i,k) in K(m,:) 
                                    jx=     t(4,2+j1,2+j2,2+j3);
                                    % add entry of the local matrix 
                                    % this row only, no duplicates if triangle
                                    if m1==i1 && m2 == i2 && m3 == i3
                                        K(m1,m2,m3,jx) = K(m1,m2,m3,jx) + Kloc(iloc,kloc);
                                        fprintf(' K(ie1%s,ie2%s,ie3%s,%2i) =   K(ie1%s,ie2%s,ie3%s,%2i) + Kloc(%2i, %2i) \n',...
                                            pm(ic1), pm(ic2), pm(ic3), jx,pm(ic1), pm(ic2), pm(ic3),jx,iloc, kloc)
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

end
function s=pm(j)
switch j
    case 1 
        s='+1';
    case 0
        s='  ';
    otherwise
        error('pm: argument must be in 0:1')
end
end