%test for vertical_sweeps

format compact
% test for vertical_sweeps
%nel=[5,4,3]
nel=[1 1 1]
n=nel+1
h = [1,1,1]
expand=1.3 
A = diag([1 1 1])
lambda=[]
params=[]
u0=[];

X = regular_mesh(nel,h,expand);
% X = add_terrain_to_mesh(X, 'hill', 'shift', 0.1)
X = add_terrain_to_mesh(X, 'hill', 'squash', 0.1)  % more thorough testing

% plot_mesh_3d(X)

% assemble sparse system matrix
[K,F,~] = sparse_assembly(A,X,u0,lambda,params);

nn = size(F,1);
%xr=rand(n1,n2,n3);
x = rand(nn,1);

K_1=nd_assembly(A,X,lambda,params);
[n1,n2,n3,m1,m2,m3]=size(K_1);
K27 = reshape(K_1,n1,n2,n3,m1*m2*m3);  % make to n x 27 nd format
K14 = ndt_convert(K27,14);

% test same results for ndt_mult from matlab and fortran
if exist('fortran/sweeps_test.exe')
    disp('testing if same result in fortran')
    err=sweeps_fortran(K,K14,n1,n2,n3,F,X,x)
    if abs(err)<1e-6
        fprintf('error %g OK\n',err)
    else
        error(sprintf('error %g too large',err))
    end
else
    warning('fortran/sweeps_test.exe not available')
end