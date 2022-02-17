function ps = make_ps(w,ts,red)
%ts is test_struct

w.tign_g = ts.c;
if ~exist('red','var')
    ps.red = subset_domain(w,1);
else
    ps.red = red;
end
new_pts = fixpoints2grid(ps.red,[ts.points(:,1),ts.points(:,2)]);

ps.paths = ts.paths;
ps.points = ts.points;
ps.grid_pts = ts.grid_pts;
ps.idx = new_pts(:,1:2);
ps.graph = ts.graph;
ps.raw_dist = ts.raw_dist;

end % function