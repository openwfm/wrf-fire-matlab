format compact
% test for ndt_f_assembly from wrfout file
path='./wrfout_d01_0001-01-01_00:00:00' % path to wrfout file
[X,u0]=read_fmw_wrfout(path);
A=diag([1 1 1]);
lambda=[];
params=[];
iflags=[2 1 1];
iflags=iflags(:);
addpath ../netcdf

% test same results for ndt_mult from matlab and fortran
% if exist('fortran/f_assembly_test.exe')
    disp('testing if same result as fortran')
    p.run_fortran=0;
    p.run_matlab=1;
    p.femwind_fortran_test=0
    err=femwind_fortran(A,X,u0,p);
    if abs(err)<1e-6
    fprintf('error %g OK\n',err)
    else
    error(sprintf('error %g too large',err))
    end
%else
%    warning('fortran/f_assembly_test.exe not available')
%end
