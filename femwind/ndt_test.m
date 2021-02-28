% test for ndt_assembly and ndt_mult
nel=[5,4,3]
nel=[1 1 1]
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
x=ones(n1,n2,n3);
y27=ndt_mult(K27,x);
K27_err_zero=big(y27)  % should be zero

disp('compare matrix vector multiply nd and ndt 27 numbers per row')
x=rand(n1,n2,n3);
y27=ndt_mult(K27,x);
yn=nd_mult(K,x);
Kn_27_err=big(yn-y27)

disp('convert to sparse compate matrix-vector multiply')
Ks = ndt_convert(K,'sparse');  
ys=Ks*x(:);
K_Ks_err=big(ys(:)-y27(:))

disp('converting to reduced storage format with 14 numbers per row')
K14 = ndt_convert(K27,14);  

disp('multiplying by all ones, should get zero')
x=ones(n1,n2,n3);
yr=ndt_mult(K14,x);
K14_err_zero=big(yr)  % should be zero
if abs(K14_err_zero)>1e-10, error('should be zero'),end
    
y14=ndt_mult(K14,x);
K14_K27_err=big(y14-y27)


K_sparse=sparse_assembly(A,X,lambda,params);
err_mat_sparse=big(Ks-K_sparse)

