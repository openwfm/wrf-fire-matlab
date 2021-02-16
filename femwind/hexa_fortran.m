function err=hexa_fortran(A,X,u0)

[Kloc,Floc,Jg]=hexa(A,X,u0)

write_array(A,'A')
write_array(X,'X')
write_array(u0,'u0')

system('./fortran/hexa_test.exe')

Kloc_f=read_array('Kloc');
Floc_f=read_array('Floc');
Jg=read_array('Jg');

Kloc_err = norm(Kloc-Kloc_f,1)
Floc_err = norm(Floc-Floc_f,1)
u0_err = norm(u0-u0_f,1)

err=max([Kloc_err,Floc_err,u0_err]);

% compare the result of matlab hexa and fortran hexa

end
