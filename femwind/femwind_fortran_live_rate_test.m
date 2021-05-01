function femwind_fortran_live_rate_test

disp('femwind_fortran_live_rate_test')
disp('run fortran testers on live data and compare')

p.fortran_test=1;
p.run_matlab=1;
femwind_rate_test(p)
end

