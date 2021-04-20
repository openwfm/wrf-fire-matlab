format compact
% test for ndt_f_assembly
nel=[5,4,3]
%nel=[1 1 1]
n=nel+1
h = [1,1,1]
expand=1.3 
A = diag([1 1 1])
%lambda=[]
params=[]
%iflags = [2 1 1]
%iflags = iflags(:)
lambda = rand(nel(1),nel(2),nel(3));

X = regular_mesh(nel,h,expand);
X = add_terrain_to_mesh(X, 'hill', 'squash', 0.1);
[CX,CH] = center_mesh(X); % get midpoints of elements
U0={rand(nel),rand(nel),rand(nel)};

% test same results for ndt_mult from matlab and fortran
if exist('fortran/ndt_w_test.exe')
disp('testing if same result in fortran')
    err=ndt_w_test_fortran(A,X, U0,lambda, params);
    if abs(err)<1e-6
    fprintf('error %g OK\n',err)
    else
    error(sprintf('error %g too large',err))
    end
else
    warning('fortran/w_assembly_test.exe not available')
end
