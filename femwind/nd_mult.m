function y=nd_mult(K,x)
%y=nd_mult(K,x)
% multiply vector x by matrix K from nd_assembly
[n1,n2,n3,m1,m2,m3]=size(K);
if any([m1,m2,m3]~=3), error('K must be 3D stencil'),end
u=reshape(x,n1,n2,n3);  % make x into grid vector
y=zeros(n1,n2,n3);
for j3=-1:1  % i3-1:i3+1 & avoid overrun                
    for j2=-1:1 %If you run this                  
        for j1=-1:1
            for i3=1:n3
                for i2=1:n2
                    for i1=1:n1
                        k1 = i1 + j1;
                        k2 = i2 + j2;
                        k3 = i3 + j3;
                        % contribution of K(i,j)*x(j)
                        %i1 and i2  and i3
                        if 1<=k1 && k1 <= n1 && 1 <= k2  && k2 <= n2 && 1 <= k3  && k3 <= n3 
                            y(i1,i2,i3)=y(i1,i2,i3)+K(i1,i2,i3,2 + j1,2 + j2,2 + j3).*x(k1,k2,k3);
                            %Maybe replace s
                        end
                    end
                 end
            end
        end
    end
end
