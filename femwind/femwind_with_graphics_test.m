function p=femwind_with_graphics_test
disp('femwind_with_graphics_test')
p=femwind_main_test
p.graphics=3; 
p.st_contour=1
p.sc_all = 1;
p.sc2_all =  [1];
p=femwind_main(p);
rate = [0.066948270621534   0.069839695192657];
if all(abs(p.rate - rate) < 1e-8)
    disp('convergence rate tests OK')
else
    error(sprintf('something changed, expected convergence rate %g %g got %g %g',...
       rate,p.rate))
end
end

