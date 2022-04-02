function [pr,pv,mpv,np] = path_ros(ps)

n = length(ps.paths);
pr = [];
pp = [];
pv = [];
for i = 2:n
    p = ps.paths(i).p;
    for j = 2:length(p)
%         d = ps.raw_dist(p(j),p(j-1));
%         t = 24*3600*(ps.points(j,3)-ps.points(j-1,3));
%         r = d/t;
%         pr = [pr;r];
        pp = [pp;p(j-1) p(j)];
 
    end
end
pp = unique(pp,'rows');
n = length(pp);
for i = 1:n

        d = ps.raw_dist(pp(i,1),pp(i,2));
        t = (ps.points(pp(i,2),3)-ps.points(pp(i,1),3));
        r = d/t/(24*3600);
        pr = [pr;r];
        rv = ros_var(d,t*24);
        pv = [pv;rv];

end
figure,histogram(pr)
pm = mean(pr);
ps = std(pr);
tstr = sprintf('ROS along paths \n Mean: %f  STD: %f',pm,ps);
title(tstr)
xlabel('ROS [m/s]')
ylabel('Frequency')

msk = isnan(pv);
nan_var = sum(msk(:));
np = nan_var/n;
mpv = mean(pv(~msk));
fprintf('Mean of std deviations on paths: %f \n',mpv)
fprintf('Percentage of path without variance: %f \n',np*100);
figure,histogram(pv(~msk));
title('Std deviation of ROS on paths')


end %function