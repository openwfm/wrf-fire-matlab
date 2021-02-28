function K2=K(K1,m2)
% K2=K(K,m2)
% convert ndt matrix from one storage format to another
[n1,n2,n3,m1]=size(K1);
t1=ndt_storage_table(m1);
t2=ndt_storage_table(m2);
K2=zeros(n1,n2,n3,m2);
for i3=1:n3
    for i2=1:n2
        for i1=1:n1
            for j3=1:3
                for j2=1:3
                    for j1=1:3
                        % location if entry (i,j) in K1  
                        k1=i1+t1(1,j1,j2,j3); % row + offset
                        k2=i2+t1(2,j1,j2,j3); % row + offset
                        k3=i3+t1(3,j1,j2,j3); % row + offset
                        kx=   t1(4,j1,j2,j3); % index in the row
                        % location if entry (i,j) in K1  
                        m1=i1+t2(1,j1,j2,j3); % row + offset
                        m2=i2+t2(2,j1,j2,j3); % row + offset
                        m3=i3+t2(3,j1,j2,j3); % row + offset
                        mx=   t2(4,j1,j2,j3); % index in the row
                        if ~(k1<1 || k1>n1 || k2<1 || k2>n2 || k3<1 || k3>n3 || ...
                             m1<1 || m1>n1 || m2<1 || m2>n2 || m3<1 || m3>n3 )
                            K2(m1,m2,m3,mx)=K1(k1,k2,k3,kx);
                        end
                    end
                end
            end
        end
    end
end

        