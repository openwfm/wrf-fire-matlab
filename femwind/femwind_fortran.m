function [W,rate]=femwind_fortran(A,X,u0,params)

W=[];

exe  = './fortran/femwind_test.exe';
write_array_nd(swap23(X{1}),'X_input');
write_array_nd(swap23(X{2}),'Y_input');
write_array_nd(swap23(X{3}),'Z_input');
write_array_nd(swap23(u0{1}),'u0_input');
write_array_nd(swap23(u0{2}),'v0_input');
write_array_nd(swap23(u0{3}),'w0_input');

write_array(A,'A_input')

% defaults
nel = size(X{1})-1;
u =zeros(nel);
v =zeros(nel);
w =zeros(nel);
rate=0;   % placeholder

disp(['running ',exe])
system(exe);
u=swap23(read_array('u'));
v=swap23(read_array('v'));
w=swap23(read_array('w'));
rate=read_array_nd('rate');
W = {u,v,w};

end


