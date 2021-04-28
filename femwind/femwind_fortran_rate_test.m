function femwind_fortran_rate_test

disp('femwind_fortran_rate_test')

p.run_fortran=1;
p.run_matlab=0;
femwind_rate_test(p)
end

