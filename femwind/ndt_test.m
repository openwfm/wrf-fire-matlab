% test for ndt_assembly and ndt_mult
nel=[5,4,3]
% nel=[1 1 1]
n=nel+1
h = [1,1,1]
expand=1.3 
A = diag([1 1 1])
lambda=[]
params=[]

X = regular_mesh(nel,h,expand);

plot_mesh_3d(X)

K=nd_assembly(A,X,lambda,params);
[n1,n2,n3,m1,m2,m3]=size(K);
K27 = reshape(K,n1,n2,n3,m1*m2*m3);  % make to n x 27 nd format

disp('testing multiplication by ones should give zero')
x1=ones(n1,n2,n3);
y27=ndt_mult(K27,x1);
K27_err_zero=big(y27)  % should be zero

disp('compare matrix vector multiply nd and ndt 27 numbers per row')
xr=rand(n1,n2,n3);
y27=ndt_mult(K27,xr);
yn=nd_mult(K,xr);
Kn_27_err=big(yn-y27)

disp('converting to reduced storage format with 14 numbers per row')
K14 = ndt_convert(K27,14); 
disp('multiplying st 14 by all ones')
x=ones(n1,n2,n3);
y1=ndt_mult(K14,x1);
K14_err_zero=big(y1)  % should be zero

disp('multiplying st 14 by random')
y14=ndt_mult(K14,xr);
K14_err_rand=big(y14-y27)  % should be zero
if abs(K14_err_rand)>1e-10, error('should be zero'),end

disp('convert to sparse compare matrix-vector multiply')
Ks = ndt_convert(K,'sparse');  
ys=Ks*xr(:);
K_Ks_err=big(ys(:)-y27(:))

K_sparse=sparse_assembly(A,X,lambda,params);
err_mat_sparse=big(Ks-K_sparse)

