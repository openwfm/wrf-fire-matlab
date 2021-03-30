function err=sweeps_fortran(K,K14,F,X,x)

% compare the result of matlab ndt_mult and fortran ndt_mult

y=vertical_sweeps(K,F,X,x);

% x_size = size(x,1);
% F_size = size(F,1);
% 
n = size(X{1});
x_r = reshape(x,n(1),n(2),n(3));
F_r = reshape(F,n(1),n(2),n(3));

write_array_nd(K14,'Kmat');
write_array_nd(F_r,'Fmat');
write_array_nd(x_r,'x_sweeps');

system('./fortran/sweeps_test.exe');

y_f=read_array_nd('x_sweeps');

err=norm(y(:)-y_f(:),inf);


end

