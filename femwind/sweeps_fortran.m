function err=sweeps_fortran(K,F,X,x)

% compare the result of matlab ndt_mult and fortran ndt_mult

y=vertical_sweeps(K,F,X,x);

write_array_nd(K,'Kmat');
write_array_nd(F,'F');
write_array_nd(x,'x_in');

system('./fortran/ndt_mult_test.exe');

y_f=read_array_nd('x_out');

err=norm(y(:)-y_f(:),inf);


end

