function err=hexa_fortran(A,X,u0)

[Kloc,Floc,Jg]=hexa(A,X,u0)

iflags = [1 1 1]
iflags = iflags(:)
write_array(A,'A')
write_array(X,'X')
write_array(u0,'u0')
write_array(iflags,'iflags')

system('./fortran/hexa_test.exe')

Kloc_f=read_array('Kloc');
Floc_f=read_array('Floc');
Jg_f=read_array('Jg');

if iflags(1) == 1
    Kloc_err = norm(Kloc-Kloc_f,1)
    err = Kloc_err;
else
    if iflags(1) == 2
        Floc_err = norm(Floc-Floc_f,1)
        err = Floc_err;
    else
        Jg_err = norm(Jg-Jg_f,1)
        err = Jg_err;
    end
end

% compare the result of matlab hexa and fortran hexa

end
