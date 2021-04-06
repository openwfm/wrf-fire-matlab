function err = ndt_assembly_err
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

%K=nd_assembly(A,X,lambda,params);
%[n1,n2,n3,m1,m2,m3]=size(K);
%K27 = reshape(K,n1,n2,n3,m1*m2*m3);

disp('converting to reduced storage format with 14 numbers per row')
%K14 = ndt_convert(K27,14); 
K = ndt_assembly(A,X,u0,lambda,params,14);

%Writing all arrays to text files for use by fortran tester
write_array_nd(A,'A');

write_array_nd(X{1},'X');
write_array_nd(X{2},'Y');
write_array_nd(X{3},'Z');


write_array_nd(K,'Kmat1');


system('./fortran/ndt_assembly_test.exe');
K_f = read_array_nd('K');

err = norm(K(:) - K_f(:),inf);
end


