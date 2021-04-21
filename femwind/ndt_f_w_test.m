format compact
addpath ../netcdf
% test for ndt_f_assembly
nel=[5,4,3]
%nel=[1 1 1]
n=nel+1
h = [1,1,1]
expand=1.3 
A = diag([1 1 1])
lambda=[]
params=[]
% [u1, u2, u3] = ndgrid()
u0={rand(nel), rand(nel), rand(nel)};
iflags = [2 1 1]
iflags = iflags(:)


X = regular_mesh(nel,h,expand);
X = add_terrain_to_mesh(X, 'hill', 'squash', 0.1);



% test same results for ndt_mult from matlab and fortran
if exist('fortran/ndt_f_test.exe')
    disp('testing if same result in fortran')
    F=ndt_f_assembly_fortran(A,X,u0,[],params);
else
    warning('fortran/ndt_f_test.exe not available')
end
