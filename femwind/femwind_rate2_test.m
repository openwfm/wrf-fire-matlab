function p=femwind_rate2_test
disp('basic convergence speed test')
p=femwind_main
p.graphics=0;
p.sc_all = 1;
p.sc2_all =  2;
p=femwind_main(p);
rate = [0 0.069839695192657];
if abs(p.rate(2) - rate(2)) < 1e-8
    disp('basic convergence rate test OK')
else
    error(sprintf('something changed, expected convergence rate %g got %g diff %g',...
       rate,p.rate,rate-p.rate))
end
end

