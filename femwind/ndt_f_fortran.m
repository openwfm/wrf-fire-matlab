function err=ndt_f_fortran(A,X,u0, iflags)

% compare the result of matlab and fortran
y= ndt_f_assembly(A,X,u0,params);

write_array(A,'A')
write_array(X,'X')
write_array(u0,'u0')
write_array(iflags,'iflags')

system('./fortran/ndt_f_test.exe');

y_f=read_array_nd('Fvec');

err= norm(y(:)-y_f(:),inf);


end
