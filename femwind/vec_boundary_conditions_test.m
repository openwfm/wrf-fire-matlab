F=rand(5,4,7);
if exist('fortran/vec_boundary_conditions_test.exe')
    disp('testing if vec_boundary_conditions gives same result in fortran')
    vec_boundary_conditions_fortran(F);
else
    warning('fortran/vec_boundary_conditions_test.exe not available')
end
