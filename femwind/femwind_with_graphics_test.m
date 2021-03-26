function femwind_with_graphics_test
disp('femwind_with_graphics_test')
p=femwind_main
p.graphics=2; 
p.st_contour=1
p.sc_all = 1;
p.sc2_all =  [1 2];
p=femwind_main(p);
rate = [0.0858361045908 0.088814805512714];
if all(abs(p.rate - rate) < 1e-8)
    disp('convergence rate tests OK')
else
    error(sprintf('something changed, expected convergence rate %g %g got %g %g',...
       rate,p.rate))
end
end

