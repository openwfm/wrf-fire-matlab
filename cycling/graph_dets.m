function graph_dets(g)

pts = [];
min_con = 9;
for i = 1:length(g)
    if sum(g(i).det(3:5)) > 0
        fires = g(i).data >= min_con;
        lons = g(i).xlon(fires);
        lats = g(i).xlat(fires);
        times = g(i).time*ones(length(lons),1);

    end
    
end


end