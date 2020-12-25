function p=hexasub(X,K,i) 
% p=hexasub(X,K,i)
% extract one submatrix from global matrix from i to i+1
% in:
%    X mesh coordinates, cell vector 3
%    K global matrix
%    i index of lower left front corner
% out:
%    p.S2 2D submatrix of K
%    p.X  3 by 8 matrix of coordinates, same as in hexa
%    p.S6 submatrix of K as 6D 
n= size(X{1});
nn = prod(n);
if any(size(K)~=nn),
    error('wrong size K')
end
p.S6=zeros(2,2,2,2,2,2);
p.S2=zeros(8,8);
p.X=zeros(3,8);
for i1=0:1
    for i2=0:1
        for i3=0:1
            ix=sub2ind(n,i(1)+i1,i(2)+i2,i(3)+i3); % global index
            ik=sub2ind([2,2,2],1+i1,1+i2,1+i3);    % local index
            p.X(1,ik)=X{1}(i(1)+i1,i(2)+i2,i(3)+i3);
            p.X(3,ik)=X{2}(i(1)+i1,i(2)+i2,i(3)+i3);
            p.X(3,ik)=X{3}(i(1)+i1,i(2)+i2,i(3)+i3);
            for j1=0:1
                for j2=0:1
                    for j3=0:1
                        jx=sub2ind(n,i(1)+j1,i(2)+j2,i(3)+j3);
                        jk=sub2ind([2,2,2],1+j1,1+j2,1+j3);
                        p.S2(ik,jk)=K(ix,jx);
                        p.S6(i1+1,i2+1,i3+1,j1+1,j2+1,j3+1)=K(ix,jx);
                    end
                end
            end
        end
    end
end
end