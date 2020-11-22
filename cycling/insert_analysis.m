function outer = insert_analysis(w,paths,a)
% puts analysis into the w file
%w = wrfout_tign(f)
%red = subset_domain(w)
%paths = graph_dets(w,1) or paths = cluster(paths,w1)
%a = squish(paths,1) , this is an anylysis fire arrival time
%                       time format is tign_g, seconds from start

red = paths.red;
if isfield(red,'red')
    red_orig = red.red;
    %interpolate back to original size grid
    F = scatteredInterpolant(red.fxlat(:),red.fxlong(:),a(:));
    a = F(red_orig.fxlat,red_orig.fxlong);
end
w.analysis = w.tign_g;
temp_a = datenum2time(a,red);
%shift possible time differneces
diff = max(w.tign_g(:))-max(temp_a);
temp_a = temp_a + diff;

w.analysis(red.ispan,red.jspan)=temp_a;

%% adapt to use new use
% restart_time=time_bounds(3);
% perimeter_time=time_bounds(4);
% wf = max(forecast - restart_time,0); % 0 in forecast fire area at restart time, >0 outside 
% wa = max(perimeter_time-analysis,0); % 0 outside of analysis fire area at perimeter time, >0 inside
% 
% % check if we have inclusion so the strip exist 
% shrink=nnz(wa + wf==0);  
% if shrink,
%     fprintf('Cannot spin up, fire area shrinking in analysis at %g points\n',shrink)
%     warning('Using analysis in place of spinup');
%     w.spinup = w.analysis
% else
% 
% % map the weights so that 0 ->1, 0->1, 
% vf=1-wf./(wf+wa);  % 1  in forecast fire area at restart time, ADDING UP TO 1 
% va=1-wa./(wf+wa); %  1  outside of analysis fire area at restart time, ADDING UP TO 1
% 
% % combine the forecast and analysis
% spinup = vf.*forecast + va.*analysis; 
% end

outer = w;
end
