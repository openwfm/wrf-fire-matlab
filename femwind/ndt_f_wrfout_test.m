format compact
% test for ndt_f_assembly from wrfout file
path = '/path/to/wrfout'
X,u0 = read_fmw_wrfout(path)
A = diag([1 1 1])
lambda=[]
params=[]
iflags = [2 1 1]
iflags = iflags(:)

% test same results for ndt_mult from matlab and fortran
if exist('fortran/ndt_f_test.exe')
    disp('testing if same result in fortran')
    err=ndt_f_fortran(A,X,u0,iflags);
    if abs(err)<1e-6
    fprintf('error %g OK\n',err)
    else
    error(sprintf('error %g too large',err))
    end
else
    warning('fortran/ndt_f_test.exe not available')
end
