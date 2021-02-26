function mt = multi_make(ps,alpha,a,b,p_param,grid_fraction)

%do make_tign on set of grids of decreasing size
gs = [2000,1000,750,600,500];
%compute grid fractions gf = gs/dlon
red = ps.red;
E = wgs84Ellipsoid;
dlon= distance(red.min_lat,red.min_lon,red.min_lat,red.max_lon,E);
dlat= distance(red.min_lat,red.min_lon,red.max_lat,red.min_lon,E);
[m,n] = size(red.tign);
gf = gs(1)/dlon;

%make series of estimates
grids = length(gs);
%10 steps, max
grids = 10;
re1 = 1;
do_plots = 0;
for i = 1:grids
    [t,e] = make_tign(ps,alpha,p_param,gf);
    if do_plots == 1
       close all
       tstr = sprintf('Iteration %d',i);
       savestr = sprintf('multi_%d',i);
       %limit plot region
       t_max = max(t(:));
       msk = t<t_max;
       min_lon = min(ps.red.fxlong(msk));
       max_lon = max(ps.red.fxlong(msk));
       lon_buff = 0.2*(max_lon-min_lon);
       min_lat = min(ps.red.fxlat(msk));
       max_lat = max(ps.red.fxlat(msk));
       lat_buff = 0.2*(max_lat-min_lat);
       figure(312)
       mesh(ps.red.fxlong(1:5:end,1:5:end),ps.red.fxlat(1:5:end,1:5:end),t(1:5:end,1:5:end));
       hold on
       scatter3(ps.points(:,2),ps.points(:,1),ps.points(:,3),'r*');
       xlim([min_lon-lon_buff max_lon+lon_buff]);
       ylim([min_lat-lat_buff max_lat+lat_buff]);
       xlabel('Lon (deg)');
       ylabel('Lat (deg)');
       zlabel('Time (datenum)');
       title(tstr);
       savefig(savestr);
       saveas(gcf,[savestr '.png']);
       hold off
    end
    ps.red.tign = t;
    mt(i).t = t;
    mt(i).e = e;
    if i < grids
        %gf = gf*gs(i)/gs(i+1);
        gf = 4/3*gf;
        close all
    end
    if i > 1
        re1 = norm(mt(i).t-mt(i-1).t)/norm(mt(i-1).t-min(mt(i-1).t(:)));
    end
    fprintf('relative error from previous estimate : \n',re1)
    if re1 < 0.01
        fprintf('not changing.....  exit loop  \n')
        break
    end
    if i == grids
        fprintf('Max steps performed\n')
    end
end


%fix the ripple on the top
r_mask = mt(end).t >=max(mt(end).t(:))-0.05;
mt(end).t(r_mask) = max(mt(end).t(:));

%compute relative error
re = norm(red.tign-mt(end).t)/(norm(red.tign-red.start_datenum));
fprintf('Relative error: %f\n',re)

end % function
