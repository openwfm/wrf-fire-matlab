function p=hexasub(X,K,i) 
% p=hexasub(X,K,i)
% extract one submatrix from global matrix from i to i+1
n= size(X{1});
p.S6=zeros(2,2,2,2,2,2);
p.S2=zeros(8,8);
for i1=0:1
    for i2=0:1
        for i3=0:1
            for j1=0:1
                for j2=0:1
                    for j3=0:1
                        ix=sub2ind(n,i(1)+i1,i(2)+i2,i(3)+i3);
                        jx=sub2ind(n,i(1)+j1,i(2)+j2,i(3)+j3);
                        ik=sub2ind([2,2,2],1+i1,1+i2,1+i3);
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