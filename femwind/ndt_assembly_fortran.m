function err=ndt_assembly_fortran(A,X,u0,lambda,params)
% call fortran version and compare results

sz = size(X{1});
disp(['ndt_assembly_fortran, mesh size ',num2str(sz)])

K_m=ndt_assembly(A,X,[],[],params,14);

%Writing all arrays to text files for use by fortran tester
write_array_nd(A,'A');
write_array_nd(swap23(X{1}),'X');
write_array_nd(swap23(X{2}),'Y');
write_array_nd(swap23(X{3}),'Z');

exe='./fortran/ndt_assembly_test.exe';
disp(['ndt_assembly done, calling ',exe])

system(exe);
K = swap23(read_array_nd('K'));

err = norm(K_m(:) - K(:),inf);

