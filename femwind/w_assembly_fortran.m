function err=w_assembly_fortran(A,X,u0,lambda,params) 
% call fortran version and compare results

lambda2 = lambda(:);
W_m=w_assembly(A,X,u0, lambda2,params);
U_m = W_m{1};
V_m = W_m{2};
W_m = W_m{3};


%Writing all arrays to text files for use by fortran tester
write_array(A,'A_test');

n = size(X{1});
lambda=reshape(lambda,n);
write_array(swap23(lambda), 'lambda_test');

write_array(swap23(X{1}),'X_test');
write_array(swap23(X{2}),'Y_test');
write_array(swap23(X{3}),'Z_test');

write_array(swap23(u0{1}),'u0_test');
write_array(swap23(u0{2}),'v0_test');
write_array(swap23(u0{3}),'w0_test');

system('./fortran/w_assembly_test.exe');
U = swap23(read_array('U_test'));
V = swap23(read_array('V_test'));
W = swap23(read_array('W_test'));



err_u = norm(U_m(:) - U(:),inf)
err_v = norm(V_m(:) - V(:),inf) 
err_w = norm(W_m(:) - W(:),inf) 
err   =  max([err_u, err_v, err_w]) 

end
