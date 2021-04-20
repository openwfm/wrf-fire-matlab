function err=ndt_assembly_fortran(A,X,u0,lambda,params); 
% call fortran version and compare resultso

lambda2 = lambda(:);
W_m=ndt_w_assembly(A,X,u0, lambda2,params);
U_m = W_m{1};
V_m = W_m{2};
W_m = W_m{3};


%Writing all arrays to text files for use by fortran tester
write_array_nd(A,'A');

write_array_nd(swap23(lambda), 'lambda');

write_array_nd(swap23(X{1}),'X');
write_array_nd(swap23(X{2}),'Y');
write_array_nd(swap23(X{3}),'Z');

write_array_nd(swap23(u0{1}),'u0');
write_array_nd(swap23(u0{2}),'v0');
write_array_nd(swap23(u0{3}),'w0');

system('./fortran/ndt_w_test.exe');
U = swap23(read_array_nd('U'));
V = swap23(read_array_nd('V'));
W = swap23(read_array_nd('W'));



err_u = norm(U_m(:) - U(:),inf)
err_v = norm(V_m(:) - V(:),inf) 
err_w = norm(W_m(:) - W(:),inf) 
err   =  max([err_u, err_v, err_w]) 

end
