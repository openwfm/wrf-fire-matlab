function K2=ndt_convert(K,m)
% K2=K(K,m2)
% convert ndt matrix from one storage format to another
switch m
    case {27,'nons'}
        m2=27;
    case {14,'symm'}
        m2=14
    case {0,'sparse'}
        m2=0;
    otherwise
        m,error('allowable values of m: 27 or nons, 14 or symm, 0 or sparse')
end
[n1,n2,n3,m1]=size(K);
t1=ndt_storage_table(m1);
if m2,
    t2=ndt_storage_table(m2);
    K2=zeros(n1,n2,n3,m);
else
    n=n1*n2*n2
    ii=zeros(n*27,1);
    jj=zeros(n*27,1);
    vv=zeros(n*27,1);
    nn=0;
end
g=@(i1,i2,i3)i1+n1*(i2-1+n3*(i3-1));  % flattened index in 1:n by 1:n matrix
for i3=1:n3
    for i2=1:n2
        for i1=1:n1
            for j3=-1:1
                for j2=-1:1
                    for j1=-1:1
                        % location of entry (i,j) in K1  
                        k1=i1+t1(1,2+j1,2+j2,2+j3); % row + offset
                        k2=i2+t1(2,2+j1,2+j2,2+j3); % row + offset
                        k3=i3+t1(3,2+j1,2+j2,2+j3); % row + offset
                        kx=   t1(4,2+j1,2+j2,2+j3); % index in the row
                        jj1 = i1+j1;
                        jj2 = i2+j1;
                        jj3 = i3+j3;
                        if  ~(jj1<1 || jj1>n1 || jj2<1 || jj2>n2 || jj3<1 || jj3>n3)
                            % j is within bounds
                            v = K(k1,k2,k3,kx); % value of the entry
                            if m2  % ndt storage scheme
                                % location of entry (i,j) in K2  
                                l1=i1+t2(1,2+j1,2+j2,2+j3); % row + offset
                                l2=i2+t2(2,2+j1,2+j2,2+j3); % row + offset
                                l3=i3+t2(3,2+j1,2+j2,2+j3); % row + offset
                                lx=   t2(4,2+j1,2+j2,2+j3); % index in the row
                                K2(l1,l2,l3,lx)=v;
                            else
                                % sparse, K2(i,j)=v
                                nn=nn+1;
                                ii(nn)=g(i1,i2,i3);
                                jj(nn)=g(jj1,jj2,jj3);
                                vv(nn)=1;
                            end
                        end
                    end
                end
            end
        end
    end
end
if mm2==0,
    K2 = sparse(ii,jj,vv,n,n);
end
end
        
