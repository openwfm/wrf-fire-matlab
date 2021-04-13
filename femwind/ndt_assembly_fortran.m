function K=ndt_assembly_fortran(A,X,u0,lambda,params,m); 
% call fortran version and compare resultso

if m~=14, 
    error('must have m=14')
end

K_m=ndt_assembly(A,X,[],[],params,m);

%Writing all arrays to text files for use by fortran tester
write_array_nd(A,'A');
write_array_nd(swap23(X{1}),'X');
write_array_nd(swap23(X{2}),'Y');
write_array_nd(swap23(X{3}),'Z');

system('./fortran/ndt_assembly_test.exe');
K = swap23(read_array_nd('K'));

err = norm(K_m(:) - K(:),inf)
end

function u=swap23(v)
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
