format compact
path='./wrfout' % path to wrfout file
out_file=sprintf('%s.mat',path);
[X,u0]=read_fmw_wrfout(path);
A=diag([1 1 1]);
lambda=[];
params=[];
iflags=[2 1 1];
iflags=iflags(:);

p=params_defaults;
p.run_fortran=0;
p.run_matlab=1;
p.femwind_fortran_test=0;
p.graphics=-2;
[W,err]=femwind_fortran(A,X,u0,p);
save(out_file)
