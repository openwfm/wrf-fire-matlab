function err = ndt_assembly_test
format compact
% test for ndt_assembly and ndt_mult
nel=[5,4,3]
%nel=[1 1 1]
n=nel+1
h = [1,1,1]
expand=1.3 
A = diag([1 1 1])
lambda=[]
params=[]
u0=[]


X = regular_mesh(nel,h,expand);
X = add_terrain_to_mesh(X, 'hill', 'squash', 0.1); 

K_m = ndt_assembly(A,X,u0,lambda,params,14);
K = ndt_assembly_fortran(A,X,u0,lambda,params,14);

err = norm(K(:) - K_m(:),inf)
end


