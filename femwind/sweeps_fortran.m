function err=sweeps_fortran(K,K14,n1,n2,n3,F,X,x)

% compare the result of matlab ndt_mult and fortran ndt_mult
y=vertical_sweeps(K,F,X,x);

% x_size = size(x,1);
% F_size = size(F,1);
% 
x_r = reshape(x,n1,n2,n3);
F_r = reshape(F,n1,n2,n3);

write_array_nd(swap23(K14),'Kmat');
write_array_nd(swap23(F_r),'Fmat');
write_array_nd(swap23(x_r),'x_sweeps');

system('./fortran/sweeps_test.exe');

y_f=read_array_nd('x_sweeps');

err=norm(y(:)-y_f(:),inf);


end

