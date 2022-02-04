function d = weak_test()
%sets up directories for weak test...

cwd = pwd;

for i = 1:9
    
    d.cores(i) = 36*i^2;
    d.nodes(i) = ceil(36*i^2/54);
    d.direct(i) = {[pwd,'/rain1_weak_ifire1_',num2str(i,'%02d')]};
    %d.direct(i) = {[pwd,'/rain1_weak_ifire2_',num2str(i,'%02d')]};
    if d.nodes(i) < 2
        queue{i} = 'small';
    else
        queue{i} = 'normal';
    end
    %     core_print = sprintf('sed -i "s/CORE_TARGET/%s/g" run_wrf_small_frontera.template',num2str(cores(i)));
    ewe = sprintf('sed -i "s/EWE/%s/g" namelist.input',num2str(42*i));
    esn = sprintf('sed -i "s/ESN/%s/g" namelist.input',num2str(42*i));
    %system call for submission script
    core_print = sprintf('sed -i "s/CORE_TARGET/%s/g" run_wrf_small_frontera.template',num2str(d.cores(i)));
    node_print = sprintf('sed -i "s/NODES_TARGET/%s/g" run_wrf_small_frontera.template',num2str(d.nodes(i)))
    %wall_print = sprintf('sed -i "s/WALL_TARGET/%s/g" run_wrf_small_frontera.template',num2str(d.wall{i}))
    queue_print = sprintf('sed -i "s/QUEUE_TARGET/%s/g" run_wrf_small_frontera.template',queue{i})
    %change queue

    
    %copy the directory
    cp_str = sprintf('cp -a ../rain1_base1 %s',d.direct{i});
    %cp_str = sprintf('cp -a ../rain1_base2 %s',d.direct{i});
    %cd to directory and change files
    system(cp_str);
    cd(d.direct{i});
    system(ewe);
    system(esn);
    
    %change submission
    system(core_print);
    system(node_print);
    system(queue_print);
    %system(wall_print);
    system('touch copied.txt');
    system('cp run_wrf_small_frontera.template run_wrf_frontera')
    
    system('./ideal.exe')
    pause(5);
    
    %cd back up
    cd(cwd)

    
end
save details.mat d


end % function