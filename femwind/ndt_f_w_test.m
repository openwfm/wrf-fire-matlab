format compact
% test for ndt_assembly and ndt_mult
nel = [4,3,2]
%nel=[1 1 1]
n=nel+1
h = [1,1,1]
expand=1.3 
A = diag([1 1 1])
lambda=[]
params=[]
u0=[];
iflags = [2 1 1]
iflags = iflags(:)

X = regular_mesh(nel,h,expand);
% X = add_terrain_to_mesh(X, 'hill', 'shift', 0.1)
X = add_terrain_to_mesh(X, 'hill', 'squash', 0.1)  % more thorough testing

% [~,F,W] = sparse_assembly(A,X,u0,lambda,params);
% F_ndt = ndt_f_assembly(A,X,params);

% F_err= norm(F(:)-F_ndt(:),inf);



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
