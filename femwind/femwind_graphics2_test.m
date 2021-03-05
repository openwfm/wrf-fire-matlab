function femwind_graphics2_test
disp('femwind_graphics2_test')
p=femwind_main
p.graphics=2; % te
p.sc_all = 1;
p.sc2_all =  [1 2];
p=femwind_main(p);
rate = [0.0858361045908 0.097188694734];
if all(abs(p.rate - rate) < 1e-8)
    disp('basic convergence rate test OK')
else
    warning(sprintf(disp('something changed, expected convergence rate %g got %g',...
       rate(2),p.rate(2)))
end
end

