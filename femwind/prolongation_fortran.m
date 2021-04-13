function u=prolongation_fortran(uc,hzc,icl3,X,params); 
% call fortran version

%Writing all arrays to text files for use by fortran tester
write_array(swap23(uc),'uc');
write_array(swap23(X{1}),'X');
write_array(swap23(X{2}),'Y');
write_array(swap23(X{3}),'Z');
write_array(hzc(:),'hcz');
write_array(icl3(:),'cl_z');

system('./fortran/prolongation_test.exe');
u = swap23(read_array('u'));

end

