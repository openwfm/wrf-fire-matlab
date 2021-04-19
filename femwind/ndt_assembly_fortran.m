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

