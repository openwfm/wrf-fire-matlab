function r=compare_fire_area(f1,f2)
% r=compare_fire_area(f1,f2)
% finds the differences in fire area between two wrfouts
% arguments:
%     f1, f2 wrfout file names
% returns:
%     structure with results

p1 = nc2struct(f1,{'Times','FIRE_AREA'},{});
p2 = nc2struct(f2,{'Times','FIRE_AREA'},{});
times1=char(p1.times)';
times2=char(p2.times)';
steps1=size(times1,1);
steps2=size(times2,1);
steps=min(steps1,steps2);
r = [];
ssum = @(a) sum(a(:));
for step=1:steps
    if(any(times1(step,:)-times2(step,:)))then
        disp(['frame ',num2str(i),' file ',f1,' time ',times1(step,:)])
        disp(['frame ',num2str(i),' file ',f2,' time ',times2(step,:)])
        warning('times in input files are not the same')
    end
    fa1 = p1.fire_area(:,:,step);  % fire areas this step
    fa2 = p2.fire_area(:,:,step);
    dfa = fa1 - fa2;
    px1 = fa1>0;                   % 0-1 arrays burning/not burning 
    px2 = fa2>0;                   % 0-1 arrays burning/not burning 
    dpx = px1 - px2;
    r.fire_area_1(1,step) = ssum(fa1);
    r.fire_area_2(1,step) = ssum(fa2);
    r.fire_area_set_diff(1,step) = ssum(abs(fa2-fa1));
    r.fire_area_diff(1,step) = ssum(fa2)--ssum(fa1);
    r.fire_area_diff_rel(1,step) = (ssum(fa2)-ssum(fa1))/((ssum(fa1)+ssum(fa1))/2+realmin);
    r.fire_area_set_diff_rel(1,step) = ssum(abs(fa2-fa1))/((ssum(fa1)+ssum(fa1))/2+realmin);
    r.fire_pixels_1(1,step) = ssum(px1);
    r.fire_pixels_2(1,step) = ssum(px2);
    r.fire_pixels_diff(1,step) = ssum(px2)-ssum(px1);
    r.fire_pixels_diff_rel(1,step) = (ssum(px2)-ssum(px1))/((ssum(px1)+ssum(px1))/2+realmin);
    r.fire_pixels_set_diff(1,step) = ssum(abs(px2-px1));
    r.fire_pixels_set_diff_rel(1,step) = ssum(abs(px2-px1))/((ssum(px1)+ssum(px1))/2+realmin);
end
end
        
        
        




