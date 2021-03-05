function K2=ndt_convert(K,st2)
% K2=K(K,st2)
% convert ndt matrix from one storage format to another
% st = storage type of output, 27 or 14 or 0=sparse
switch st2
    case {27,'nons'}
        st2=27;
    case {14,'symm'}
        st2=14
    case {0,'sparse'}
        st2=0;        
    otherwise
        st2,error('allowable values of st: 27 or nons, 14 or symm, 0 or sparse')
end
[n1,n2,n3,st1]=size(K);
t1=ndt_storage_table(st1);
if st2, % ndt storage
    t2=ndt_storage_table(st2);
    K2=zeros(n1,n2,n3,st2);
else  % sparse storage
    n=n1*n2*n3;
    mnz=n*27; % max nonzeros
    ii=zeros(mnz,1);
    jj=zeros(mnz,1);
    vv=zeros(mnz,1);
    nn=0;
    % K2=sparse([],[],[],n,n,mnz);
end
g=@(i1,i2,i3)i1+n1*(i2-1+n2*(i3-1));  % flattened index in 1:n by 1:n matrix
for i3=1:n3
    for i2=1:n2
        for i1=1:n1
            for j3=-1:1
                for j2=-1:1
                    for j1=-1:1
                        % global index of neighbor node k
                        k3 = i3+j3;
                        k2 = i2+j2;
                        k1 = i1+j1;
                        % relative index of the row where the entry (i,k) 
                        % is stored, 0 0 0 = this row 
                        % in fortran we won't have the 2+ because
                        % the array t will be indexed -1:1
                        % row m where the entry (i,k) is stored in K
                        m3 = i3+t1(3,2+j1,2+j2,2+j3);
                        m2 = i2+t1(2,2+j1,2+j2,2+j3);
                        m1 = i1+t1(1,2+j1,2+j2,2+j3);
                        % index of the matrix entry (i,k) in K(m,:) 
                        jx = t1(4,2+j1,2+j2,2+j3);
                        if ~(m1<1 || m1>n1 || m2<1 || m2>n2 || m3<1 || m3>n3 || ...
                             k1<1 || k1>n1 || k2<1 || k2>n2 || k3<1 || k3>n3 )   
                           % j is within bounds
                            v = K(m1,m2,m3,jx); % value of the entry
                            if st2  % ndt storage scheme
                                % entry (i,k) is stored in K2(l,p)  
                                l3=i3+t2(3,2+j1,2+j2,2+j3); % row + offset
                                l2=i2+t2(2,2+j1,2+j2,2+j3); % row + offset
                                l1=i1+t2(1,2+j1,2+j2,2+j3); % row + offset
                                if ~(l1<1 || l1>n1 || l2<1 || l2>n2 || l3<1 || l3>n3) 
                                    px = t2(4,2+j1,2+j2,2+j3); % index in the row
                                    K2(l1,l2,l3,px)=v;
                                end
                            else % sparse
                                i=g(i1,i2,i3);
                                k=g(k1,k2,k3);
                                % K2(i,k)=v;
                                nn=nn+1;
                                ii(nn)=i;
                                jj(nn)=k;
                                vv(nn)=v;
                            end
                        end
                    end
                end
            end
        end
    end
end
if st2==0
    K2 = sparse(ii(1:nn),jj(1:nn),vv(1:nn),n,n);
end
end
        
