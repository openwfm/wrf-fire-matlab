format compact
% test for ndt_assembly and ndt_mult
nel = [5,4,3]
%nel=[1 1 1]
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

plot_mesh_3d(X)

K=nd_assembly(A,X,lambda,params);
[n1,n2,n3,m1,m2,m3]=size(K);
K27 = reshape(K,n1,n2,n3,m1*m2*m3);  % make to n x 27 nd format

disp('converting to reduced storage format with 14 numbers per row')
K14 = ndt_convert(K27,14); 
disp('multiplying st 14 by all ones')
x1=ones(n1,n2,n3);
y1=ndt_mult(K14,x1,3);
K14_err_zero=big(y1)  % should be zero
if abs(K14_err_zero)>1e-10, error('should be zero'),end

disp('testing multiplication by ones should give zero')
x1=ones(n1,n2,n3);
y27=ndt_mult(K27,x1);
K27_err_zero=big(y27)  % should be zero


disp('compare matrix vector multiply nd and ndt 27 numbers per row')
xr=rand(n1,n2,n3);
y27=ndt_mult(K27,xr);
yn=nd_mult(K,xr);
Kn_27_err=big(yn-y27)

y14=ndt_mult(K14,xr);
K14_err_rand=big(y14-y27)  % should be zero
if abs(K14_err_rand)>1e-10, error('should be zero'),end

disp('convert to sparse compare matrix-vector multiply')
Ks = ndt_convert(K,'sparse');  
ys=Ks*xr(:);
K_Ks_rand_err=big(ys(:)-y27(:))
if abs(K_Ks_rand_err)>1e-10, error('should be zero'),end

if exist('fortran/ndt_boundary_conditions_test.exe')
    disp('testing if ndt_boundary_conditions gives same result in fortran')
    ndt_boundary_conditions_fortran(K14);
else
    warning('fortran/ndt_boundary_conditions_test.exe not available')
end

disp('testing multiply vec by st 14 vs st 27')

% test same results for ndt_mult from matlab and fortran
if exist('fortran/ndt_mult_test.exe')
    disp('testing if ndt_mult same result in fortran')
    err=ndt_mult_fortran(K14,xr);
    if abs(err)<1e-6
        fprintf('error %g OK\n',err)
    else
        error(sprintf('error %g too large',err))
    end
else
    warning('fortran/ndt_mult_test.exe not available')
end

disp('sparse assembly vs ndt_assembly test')
K_sparse=sparse_assembly(A,X,u0,lambda,params);
err_mat_sparse=big(Ks-K_sparse);
if abs(err_mat_sparse)>1e-10, error('should be zero'),end

K27a = ndt_assembly(A,X,u0,lambda,params,27);
K27a_err = big(K27 - K27a)
if abs(K27a_err)>1e-10, error('should be zero'),end


K14a = ndt_assembly(A,X,u0,lambda,params,14);
disp('multiplying ndt_assembly st 14 by all ones')
x=ones(n1,n2,n3);
y1=ndt_mult(K14,x1);
K14a_err_zero=big(y1)  % should be zero
if K14a_err_zero > 1e-6
    error('residual should be zero')
end

disp('comparing st 14 with converged st 27')
K14a_err = big(K14 - K14a)
if abs(K14a_err)>1e-10, error('should be zero'),end

disp('ndt_assembly_fortran')
err = ndt_assembly_fortran(A,X,u0,lambda,params,14)
if err > 1e-6
    error('ndt_assembly_fortran error too large')
end
