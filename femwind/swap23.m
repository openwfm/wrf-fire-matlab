function u=swap23(v)
% swap indices to/from WRF ikj ordering
% arrays of rank up to 4
[n1,n2,n3,n4]=size(v);
u=zeros(n1,n3,n2,n4);
for l=1:n4
    for k=1:n3
        for j=1:n2
            for i=1:n1
                 u(i,k,j,l)=v(i,j,k,l);
            end
        end
    end
end
end

