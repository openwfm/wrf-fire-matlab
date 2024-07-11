function femwind_run_fortran_rate_test

disp('femwind_fortran_rate_test')
disp('run fortran version only')

p.nelem3=[2500,2500,10];
p.run_fortran=1;
p.run_matlab=0;
femwind_rate_test(p)
end

