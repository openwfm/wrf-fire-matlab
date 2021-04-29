format compact
% test for ndt_f_assembly from wrfout file
path='./wrfout' % path to wrfout file
[X,u0]=read_fmw_wrfout(path);
A=diag([1 1 1]);
lambda=[];
params=[];
iflags=[2 1 1];
iflags=iflags(:);

p=params_defaults;
p.run_fortran=0;
p.run_matlab=1;
p.femwind_fortran_test=0
[W,err]=femwind_fortran(A,X,u0,p);
if abs(err)<1e-6
    fprintf('error %g OK\n',err)
else
    error(sprintf('error %g too large',err))
end
