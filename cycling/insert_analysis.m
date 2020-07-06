function outer = insert_analysis(w,paths,a)

%w = wrfout_tign(f)
%red = subset_domain(w)
%paths = graph_dets(w,1)
%a = squish(paths,1)

red = paths.red;
w.analysis = w.tign_g;
temp_a = datenum2time(a,red);
diff = max(w.tign_g(:))-max(temp_a);p2,
temp_a = temp_a + diff;
w.analysis(red.ispan,red.jspan)=temp_a;

outer = w;
end