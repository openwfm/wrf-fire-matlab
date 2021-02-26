function y=nd_mult(K,x)
%y=nd_mult(K,x)
% multiply vector x by matrix K from nd_assembly
[n1,n2,n3,m1,m2,m3]=size(K);
if any([m1,m2,m3]~=3), error('K must be 3D stencil'),end
u=reshape(x,n1,n2,n3);  % make x into grid vector
y=zeros(n1,n2,n3);
for i3=1:n3
    for i2=1:n2
        for i1=1:n1
            s=0;
            for j3=max(i3-1,1):min(i3+1,3)  % i3-1:i3+1 & avoid overrun
                for j2=max(i2-1,1):min(i2+1,3)
                    for j1=max(i1-1,1):min(i1+1,3)
                        % contribution of K(i,j)*x(j)
                        i1,i2,i3,j1,j2,j3
                        % print('i1 =', i1, 'i2 =',i2,'i3 =',i3,'j1 =',j1,'j2 =',j2,'j3 =',j3)  
                        s=s+K(i1,i2,i3,2+j1-i1,2+j2-i2,2+j3-i3)*x(j1,j2,j3);
                    end
                end
            end
            y(i1,i2,i3)=s;
        end
    end
end
