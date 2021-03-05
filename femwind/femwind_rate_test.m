function femwind_rate_test
disp('basic convergence speed test')
p=femwind_main
p.sc_all = 1
p.sc2_all = 1
p=femwind_main(p)
rate = 0.0858361045908
if abs(p.rate - rate) < 1e-8
    disp('basic convergence rate test OK')
else
    disp('something changed, expected convergence rate %g got %g',...
       rate,p.rate)
end
end

