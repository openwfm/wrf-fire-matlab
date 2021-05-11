function err=ndt_mult_fortran(kmat,u)

% compare the result of matlab ndt_mult and fortran ndt_mult

y=ndt_mult(kmat,u);

write_array_nd(kmat,'kmat');
write_array_nd(u,'u');

system('./fortran/ndt_mult_test.exe');

y_f=read_array_nd('y');

err=norm(y(:)-y_f(:),inf)


end
