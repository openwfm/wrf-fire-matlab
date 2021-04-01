function femwind_bench_test
disp(' convergence speed test')
p=femwind_main
p.graphics=0;
p.sc_all = [1 2];
p.sc2_all =  [2 4 8 16 32];
p.save_files=3;
p.levels=6
p.save_file_prefix='bench_6levels';
p=femwind_main(p);
end

