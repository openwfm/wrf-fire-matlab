% test for nd_assembly and nd_mult
nel=[5,4,3]
h = [1,1,1]
expand=1 
A = diag([1 1 1])
lambda=[]
params=[]

X = regular_mesh(nel,h,expand);

K=nd_assembly(A,X,lambda,params);

x=ones(size(X{1}));
!git 
y=nd_mult(K,x);

err_zero=big(y)  % should be zero