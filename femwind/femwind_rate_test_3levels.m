function p=femwind_rate_test(pin)

disp('basic convergence speed test')
p=femwind_main;
p.graphics=-2;
p.sc_all = 1;
p.sc2_all = 1;
% p.levels=8;
% params.P_by_x=1;  % coarsening proportioned by x  % now always
p.coarse_K=1; % 'variational';  
p.coarse_K=2; % 'assembly';
p.levels=15;
p.levels=3; fprintf('levels=%i for debugging only!',p.levels)
p.debug_level=0;
p.test_fortran=0;
if exist('pin','var')
    f=fieldnames(pin);;
else
    f={};
end
for i=1:length(f)
    ff=f{i};
    fprintf('overiding p.%s=',ff)
    disp(pin.(ff))
    p.(ff)=pin.(ff);
end
p
p=femwind_main(p);
rates(3,1) =  0.066948270621534;  % 3 levels, coarse P variational
rates(8,1) =  0.066931756926231;  % 8 levels, coarse P variational
rates(8,2) =  0.066930988019406;  % 8 levels, coarse P assembly
rates(15,2)=  0.066930988013484;  % max levels, coarse P assembly
rate = rates(p.levels,p.coarse_K);
if abs(p.rate - rate) < 1e-8
    disp('basic convergence rate test OK')
else
    error(sprintf('something changed, expected convergence rate %g got %g diff %g',...
       rate,p.rate,rate-p.rate))
end
end

