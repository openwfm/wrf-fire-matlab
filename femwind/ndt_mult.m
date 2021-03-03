function y=ndt_mult(K,x)
%y=nd_mult(K,x)
% multiply vector x by matrix K from nd_assembly
[n1,n2,n3,m]=size(K);
t = ndt_storage_table(m); 
u=reshape(x,n1,n2,n3);  % make x into grid vector if needed
y=zeros(n1,n2,n3);
n=0;
for i3=1:n3   % global index of row = node i
    for i2=1:n2
        for i1=1:n1
            for j3=-1:1   % relative index j of neighbor = nonzero in the row i               
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
                        % row m where the entry (i,k) is stored
                        m3 = i3+t(1,2+j1,2+j2,2+j3);
                        m2 = i2+t(1,2+j1,2+j2,2+j3);
                        m1 = i1+t(1,2+j1,2+j2,2+j3);
                        % index of the matrix entry (i,k) in K(m,:) 
                        jx=t(4,2+j1,2+j2,2+j3);
                        % contribution of K(i,j)*x(j) if index not out of bounds
                        % in fortran, we won't need to worry about out of bounds
                        % because K and x will be wrapped with zeros
                        if ~(m1<1 || m1>n1 || m2<1 || m2>n2 || m3<1 || m3>n3 || ...
                             k1<1 || k1>n1 || k2<1 || k2>n2 || k3<1 || k3>n3 ) 
                            n=n+1;
                            v=K(m1,m2,m3,jx);s=u(k1,k2,k3);
%                            fprintf('i=(%g %g %g) j=(%g %g %g) m=(%g %g %g) jx=%g k=(%g %g %g) v=%g u=%g\n',...
%                                     [i1 i2 i3     j1 j2 j3     m1 m2 m3     jx    k1 k2 k3     v    s])
                            fprintf('i=(%g %g %g) j=(%g %g %g) m=(%g %g %g) k=(%g %g %g) v=%g u=%g\n',...
                                     [i1 i2 i3     j1 j2 j3     m1 m2 m3       k1 k2 k3     v    s])
                            y(i1,i2,i3)=y(i1,i2,i3)+v*s;
                        end
                    end
                 end
            end
        end
    end
end
