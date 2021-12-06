disp('before merging all tests here must run successfully')
nd_test
ndt_test
hexa_test
nd_test
ndt_assembly_test
f_assembly_test
w_assembly_test
vec_boundary_conditions_test
vertical_test
p=[];
p.test_fortran=1;
femwind_rate_test(p)
