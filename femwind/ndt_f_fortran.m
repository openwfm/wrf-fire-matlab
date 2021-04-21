function err=ndt_f_fortran(A,X,u0,iflags)

% compare the result of matlab and fortran

write_array_nd(swap23(X{1}),'X');
write_array_nd(swap23(X{2}),'Y');
write_array_nd(swap23(X{3}),'Z');
write_array_nd(swap23(u0{1}),'Xu0');
write_array_nd(swap23(u0{2}),'Yu0');
write_array_nd(swap23(u0{3}),'Zu0');

write_array(A,'A')
write_array(iflags,'iflags')

disp('running ./fortran/ndt_f_test.exe')
system('./fortran/ndt_f_test.exe');

y_f=swap23(read_array_nd('Fvec'));

disp('running ndt_f_assembly.m')

y= ndt_f_assembly(A,X,u0);

err= norm(y(:)-y_f(:),inf);

end
