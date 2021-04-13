function uc=restriction_fortran(u,hzc,icl3,X,params); 
% call fortran version

% get coarse size
n=size(X{1});
[icl1,icl2]=hzc2icl(hzc,n);
nc = [length(icl1),length(icl2),length(icl3)];

%Writing all arrays to text files for use by fortran tester
uc = rand(nc);  ! for sizes only
write_array(swap23(uc),'uc');
write_array(swap23(u),'u');
write_array(swap23(X{1}),'X');
write_array(swap23(X{2}),'Y');
write_array(swap23(X{3}),'Z');
write_array(hzc(:),'hcz');
write_array(icl3(:),'cl_z');

system('./fortran/restriction_test.exe');
uc = swap23(read_array('uc'));

end

