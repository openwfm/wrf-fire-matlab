function y=nd_mult(K,x)
%y=nd_mult(K,x)
% multiply vector x by matrix K from nd_assembly
[n1,n2,n3,m1,m2,m3]=size(K);
if any([m1,m2,m3]~=3), error('K must be 3D stencil'),end
u=reshape(x,n1,n2,n3);  % make x into grid vector
y=zeros(n1,n2,n3);
for j3=-1:1                  
    for j2=-1:1               
        for j1=-1:1
            for i3=max(1,1-j3):min(n3, n3-j3) 
                k3 = i3 + j3;                             
                for i2=max(1,1-j2):min(n2, n2-j2)
                    k2 = i2 + j2;
                    for i1=max(1,1-j1):min(n1, n1-j1)
                        k1 = i1 + j1;
                        % contribution of K(i,j)*x(j)
                        %i1 and i2  and i3
                        y(i1,i2,i3)=y(i1,i2,i3)+K(i1,i2,i3,2 + j1,2 + j2,2 + j3)*x(k1,k2,k3);                       
                    end
                 end
            end
        end
    end
end
