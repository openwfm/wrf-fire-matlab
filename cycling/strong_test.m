function details = strong_test()
%sets up directories for strong scaling test

for i = 1:20
    
    cores(i) = (2*i -1)*54;
    nodes(i) = 2*i - 1;
    wall{i} = '24:00:00'
    wall_time(i) = 24;
    if nodes(i) > 1
       wall_time(i) =  ceil(24/(nodes(i)*.5));
       wall{i} = sprintf('%02d:00:00',wall_time(i))
    end
    %change queue
    if nodes < 2
        queue{i} = 'small';
    else
        queue{i} = 'normal';
    end
    
    direct(i) = {[pwd,'/ff_sfc_strong_',num2str(i,'%04d')]};
    
    %mkdir_print = sprintf('mkdir %s',direct{i});
    core_print = sprintf('sed -i "s/CORE_TARGET/%s/g" run_wrf_small_frontera.template',num2str(cores(i)));
    node_print = sprintf('sed -i "s/NODES_TARGET/%s/g" run_wrf_small_frontera.template',num2str(nodes(i)))
    wall_print = sprintf('sed -i "s/WALL_TARGET/%s/g" run_wrf_small_frontera.template',num2str(wall{i}))
    queue_print = sprintf('sed -i "s/QUEUE_TARGET/%s/g" run_wrf_small_frontera.template',queue{i})
    cp_print = sprintf('cp -a ../fireflux_sfc_strong_base %s',direct{i});
    system(cp_print);
    
    
    %change submission script
    cd(direct{i});
    system(core_print);
    system(node_print);
    system(queue_print);
    system(wall_print);
    system('cp run_wrf_small_frontera.template run_wrf_small_frontera')
    
    save details.mat 
    
    %move back up
    cd ..
    
    
    
    
    
end %if i = ...

[nodes',cores']

details.direct = direct';
details.cores = cores';
details.nodes = nodes';
details.wall = wall';
details.wall_time = wall_time';

save details.mat details

end % function