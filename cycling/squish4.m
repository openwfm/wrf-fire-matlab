function tign_new = squish4(ps,a,b,gq)
%% ps is struct with paths, red, graph,distances, etc.
% a, b - parameters of the function smooth up function
%% ps = graph_dets(w)

if ~exist('a','var')
    a = 100;
    b = 1/3;
end

%combine l2 time with detection points fixed to grid
pts = ps.grid_pts;
pts(:,3) = ps.points(:,3);
pts(1,3) = pts(1,3)-1/2;
%forecast
tign = ps.red.tign;
% t_max = max(tign(:))-.1;
% tign(tign>=t_max) =t_max+1;
t0 = min(tign(:));
[m,n]=size(ps.red.tign);

% flat_tign = input_num('Use flat start? 1 = yes',0)
% if flat_tign
%     tign = ps.red.end_datenum*ones(m,n);
% end

%%% data likelihood spline %%%
%set up data likelihood spline
[p_like_spline,~,~] = make_spline(100,1000);
%test
% t = 4*(-10:20)*3600;
% ts = p_like_spline(t);
% figure,plot(t,exp(ts))
% hold on,plot(t,1-exp(ts))
%make vector of data likelihoods
if ~exist('gq','var')
    gq = input_num('use ground detections? yes = [1]',1,1);
end
if gq
    smooth_ground = 1;%max(m,n)/50;
    fprintf('Collecting ground detection data\n')
    %grounds = ground_detects(ps.red);
    %[jgrid,igrid]=meshgrid([1:length(ps.red.jspan)]',[1:length(ps.red.ispan)]');
    [jgrid,igrid]=meshgrid([1:n]',[1:m]');
    %make polygon around detections
    %infire = inpolygon(grounds.land(:,4),grounds.land(:,3),pts(:,2),pts(:,1));
    %infire = inpolygon(grounds.land(:,4),grounds.land(:,3),grounds.land(infire,4),grounds.land(infire,3));
    infire = inpolygon(igrid(:),jgrid(:),ps.idx(:,1),ps.idx(:,2));
    %remove holes in mask
    infire = inpolygon(igrid(:),jgrid(:),igrid(infire),jgrid(infire));
end
%try

tmax=max(tign(:));
tign_ground = tign;
tign_flat = ps.red.end_datenum*ones(size(ps.red.tign));
%figure(159),hold on,scatter(pts(:,2),pts(:,1),'*r')
time_err = 0.0;
g_diff = (tign_flat-tign+time_err)*24*3600;
g_likes = p_like_spline(g_diff);
beta = 1/10;
beta_vect = exp(g_likes);
use_beta_likes = 1;
%%% ground detection likelihood
% g_times = zeros(pts_length,1);
% for i = 1:pts_length
%     g_times(i) = ps.red.tign(ps.idx(i,1),ps.idx(i,2));
% end
ground_steps = 1;
if gq
    data_area = sum(infire);
    for i = 1:ground_steps
        fire_mask = tign_ground < tmax-0.1;
        fire_area = sum(fire_mask(:));
        out_side = fire_mask(:)-infire;
        out_side(out_side<0)=0;
        out_sum(i) = sum(out_side(:));
        out_fraction = out_sum(i)/data_area;
        if i > 2 && out_sum(i) == out_sum(i-1)
            smooth_ground = 0.9*smooth_ground;
            fprintf('decreasing ground smoothing \n')
        end
        if out_fraction < 0.1
            break
        end
        fprintf('Fire area: %f pixel_out/data area: %f  pixels_out %f\n',fire_area,out_sum(i)/data_area,out_sum(i));
        if use_beta_likes == 1
            tign_ground(~infire) = beta_vect(~infire).*tign_ground(~infire)+(1-beta_vect(~infire)).*tign_flat(~infire);
        else
            tign_ground(~infire) = beta*tign_ground(~infire)+(1-beta)*tign_flat(~infire);
        end
        %flatten the top
        t_mask = tign_ground > tmax;
        tign_ground(t_mask) = tmax;
        %tign_temp = imgaussfilt(tign_ground,smooth_ground );
        tign_temp = smooth_up(ps.red.fxlong,ps.red.fxlat,tign_ground,a,b);
        tign_ground(~infire) = tign_temp(~infire);
%         a = 1/10;%1/2-1/(2*i);
%         tign_ground(infire) = a*tign_ground(infire)+(1-a)*tign_temp(infire);
%         figure(73),scatter3(ps.red.fxlong(~infire),ps.red.fxlat(~infire),tign_ground(~infire))
%         pause(1/2)
        figure(159)
        contourf(ps.red.fxlong,ps.red.fxlat,tign_ground,20,'k'),hold on
        scatter(pts(1:10:end,2),pts(1:10:end,1),'*r'),hold off
        t_str = sprintf('Perimeter Shrinking \n Iteration %d',i);
        save_str = sprintf('perim_shrink_%d',i);
        xlabel('Lon'),ylabel('Lat'),title(t_str)
%         xl=[-121.5094 -121.2496];xlim(xl)
%         yl = [46.0367   46.2035];ylim(yl)
%         savefig(save_str);
%         saveas(gcf,[save_str '.png']);
        %figure(160),mesh(ps.red.fxlong,ps.red.fxlat,tign_ground);
        %pause(.5)
        %t_min(i) = min(tign_ground(:));
    end
end
%plot(t_min)
% figure,scatter(igrid(infire),jgrid(infire))
% time_mask = grounds.land(:,5) > max(ps.points(:,3));
% ground_mask = logical(time_mask.*~infire);
tign_new = tign_ground;
%%%plot and save original perimeter
% figure(159),hold on
% scatter(pts(:,2),pts(:,1),'*r')
% contourf(ps.red.fxlong,ps.red.fxlat,ps.red.tign,20,'k'),hold off
% t_str = sprintf('Perimeter Shrinking \n Iteration %d',0);
% save_str = sprintf('perim_shrink_%d',0);
% xlabel('Lon'),ylabel('Lat'),title(t_str)
% xlim(xl);
% ylim(yl);
% savefig(save_str);
% saveas(gcf,[save_str '.png']);

%%%%%%% end ground detection work %%%%%%


idx = ps.idx;
fig_num = 33;
pts_length = length(ps.grid_pts);
%max time to look at detection data
max_l2_time = max(max(pts(:,3)));




%tign_flat is constant fire arrival time
%tign_flat=ps.red.end_datenum*ones(size(tign));
title_str = 'Analysis';

%make matrix of forecast tign time for the points in the graph
t_times = zeros(pts_length,1);
for i = 1:pts_length
    t_times(i) = tign(idx(i,1),idx(i,2));
end



%data likelikehoods
time_diff = (pts(:,3)-t_times)*24*3600;
likes = p_like_spline(time_diff);
alpha_vect = exp(likes);

%mask = ones(size(tign));
%change back to z = tign
% if exist('ps.tn','var')
%     figure(68),contour3(ps.red.fxlong,ps.red.fxlat,ps.tn-ps.red.start_datenum,20)
% else
%     figure(68),contour3(ps.red.fxlong,ps.red.fxlat,tign-ps.red.start_datenum,20)
% end
% hold on
%lons = ps.red.fxlong;
%lats = ps.red.fxlat;

%just squish the ignition point
%bigger k ==> more smoothing

norms=[];

%random multiplier, increase for larger grids
%perturbs points on path in x-y plane
% try computing this as a fraction of grid size
rm = round(m/100)+1;

% random multiplier, keep the same
% perturbs points downward in time to
rt = 1;
% weight for tign_new

%alhpa blends  estimate of tign at a point with old estimate
% new_estimate = alph*current_setimate + (1-alpha)*old_estimate
% for data assimilation, alpha will be computed from exp(likelihood)
% alpha = 0.5;
%constant for smooth in rlx_shp
alpha_2 = 0.7; %smaller alph_2 ==> smoother
%number of loops to run
smoothings = 20;
for k = 1:smoothings
%     figure(fig_num),mesh(ps.red.fxlong,ps.red.fxlat,tign_new)
%     title(title_str)
    %hold on,scatter3(ps.grid_pts(:,2),ps.grid_pts(:,1),ps.points(:,3),'*r'),hold off
    %pause(3/k)
    %plotting a slice
    %figure(fig_num+17),plot(tign_new(:,202));
    for i = 1:length(ps.paths)
        p = ps.paths(i).p;
        %         figure(73),hold on
        %         plot3(pts(p,2),pts(p,1),pts(p,3)-ps.red.start_datenum,'r')
        %         hold off
        %plot3(pts(p,2),pts(p,1),tign(idx(p,1),idx(p,2))-ps.red.start_datenum,'g')
        for j = 1:length(p) %p(j) is the current detection vertex
            tign_old = tign_new;
            %mesh indices for path points, perturbed
            p_i = idx(p(j),1);%+rm*round(randn);
            p_j = idx(p(j),2);%+rm*round(randn);
            %%% make mean of old and new, in small block around path point
            %alpha is now the data likelikehood
            %alpha = alpha_vect(p(j));
            alpha = 0.0;
            tign_new(p_i-round(rm*rand):p_i+round(rm*rand),p_j-round(rm*rand):p_j+round(rm*rand)) = alpha*tign_new(p_i,p_j) + (1-alpha)*pts(p(j),3)-rt;
            %%%% alternate strategy.
            %             if k == 1 && j > 1
            %             t1 = min(tign_new(p_i,p_j),tign_new(p_i-1,p_j-1));
            %             t2 = max(tign_new(p_i,p_j),tign_new(p_i-1,p_j-1));
            %             dt12 = t2-t1;
            %             dg = pts(p(j),3)-pts(p(j-1),3);
            %             time_shift = 0.5*(dg-dt12);
            %             t1 = t1-time_shift;
            %             t2 = t2+time_shift;
            %             tign_new(p_i,p_j) = max(t1,t2);
            %             tign_new(p_i-1,p_j-1) = min(t1,t2);
            %             end
            
            %tign_flat(p_i-round(rm*rand):p_i+round(rm*rand),p_j-round(rm*rand):p_j+round(rm*rand)) = 0.5*(tign_new(p_i,p_j) + pts(p(j),3)-rt*rand);
            
            %no need to interpolate - already done
            % interpolate new point between adjacent path detections
            if j > 1
                %weighted average will move new point close to that with
                %higher FRP
                %                 frp1 = ps.points(p(j-1),5);
                %                 frp2 = ps.points(p(j),5);
                %                 w1 = frp1/(frp1+frp2);
                %                 w2 = frp2/(frp1+frp2);
                w1 = 0.5;
                w2 = 0.5;
                %assign tign for all in small, block around midpoint
                %             frp1 = ps.points(p(j-1),5);
                %             frp2 = ps.points(p(j),5);
                %             w1 = frp1/(frp1+frp2);
                %             w2 = frp2/(frp1+frp2);
                new_lat = w1*pts(p(j-1),1)+w2*pts(p(j),1);
                new_lon = w1*pts(p(j-1),2)+w2*pts(p(j),2);
                new_t = w1*pts(p(j-1),3)+w2*pts(p(j),3);
                [new_i,new_j,new_lat,new_lon]= fixpt(ps.red,[new_lat,new_lon]);
                tign_new(new_i-round(rm*rand):new_i+round(rm*rand),new_j-round(rm*rand):new_j+round(rm*rand)) = new_t-rt*rand;
                
                %tign_new(new_i-round(rm*rand):new_i+round(rm*rand),new_j-round(rm*rand):new_j+round(rm*rand)) = new_t-rt*rand;
                
                %mask(new_i,new_j)=0;
            end
        end
    end
    %size of local averaging to apply aoutomate by grid size?
    %patch = max(1,round(sqrt(smoothings-k)));
    patch = 4*ceil(max(m,n)/100);
    fprintf('Using patch size %d in rlx_shp function \n',patch);
    %patch = 2;
    
    %%%%% make use of ground detections %%%%
    %tign_new(~infire) = beta*tign_new(~infire)+(1-beta)*tign_flat(~infire);
    
    %smooth the tign
    %tign_new(tign_new < t0) = t0;
    tign_new = smooth_up(ps.red.fxlong,ps.red.fxlat,tign_new,a,b);
    if k < smoothings/5
        tign_new = imgaussfilt(tign_new,1);
        tign_new = rlx_shp(tign_new,alpha_2,patch);
    end
%   tign_flat = rlx_shp(tign_flat,alpha_2,patch);
    
    %collect information about tign at the detection points
    for i = 1:pts_length
        t_times(i) = tign_new(ps.idx(i,1),ps.idx(i,2));
        %flat_times(i) = tign_flat(ps.idx(i,1),ps.idx(i,2));
    end
    norms(k,1) = norm(t_times-ps.points(:,3),2);
    norm_tign_new = norms(k,1);

    
    %only do norm for times before final detection time
    time_mask = tign_new < pts(end,3);  %max(max(pts(:,3)));
    norms(k,2) = norm(tign_new(time_mask)-tign_old(time_mask));%/norm(tign_old(time_mask));
    figure(fig_num+3);plot(norms(1:k,1));
    tstr = sprintf('Norm of difference between successive \n TIGN after each interpolation');
    title(tstr)
    figure(fig_num+4);plot(1:k,norms(1:k,2))
    tstr = sprintf('Norm of difference between times of \n detections and TIGN at detectionon locations');
    title(tstr);
    xlabel('Iterations of Interpolation')
    %temp_var(k) = min(tign_new(:))-ps.red.start_datenum;
    %figure(fig_num+5);plot(1:k,temp_var(1:k)*24),title('Change in ignition time'),xlabel('iteration'),ylabel('hours')
    fprintf('Loop %d complete norm of diff = %f \n', k,norms(k))
    norm_break = 1;
    if k > 2 && norms(k,norm_break) >= norm(k-1,norm_break) && norms(k,norm_break) >= norms(k-2,norm_break)
        fprintf('graph norm increase \n')
        rm = 0;
        rt = 0;
        %break
    end
end
%tign_new = [];
%smooth the result


% %%% try using ground detections after paths....
ground_steps = 4;
tign_ground = tign_new;
if gq
    data_area = sum(infire);
    for i = 1:ground_steps
        fire_mask = tign_ground < tmax-0.1;
        fire_area = sum(fire_mask(:));
        out_side = fire_mask(:)-infire;
        out_side(out_side<0)=0;
        out_sum(i) = sum(out_side(:));
        out_fraction = out_sum(i)/data_area;
        if i > 2 && out_sum(i) == out_sum(i-1)
            smooth_ground = 0.9*smooth_ground;
            fprintf('decreasing ground smoothing \n')
        end
        if out_fraction < 0.1
            %break
        end
        fprintf('Fire area: %f pixel_out/data area: %f  pixels_out %f\n',fire_area,out_sum(i)/data_area,out_sum(i));
        tign_ground(~infire) = beta*tign_ground(~infire)+(1-beta)*tign_flat(~infire);
        %tign_ground(~infire) = beta_vect(~infire).*tign_ground(~infire)+(1-beta_vect(~infire)).*tign_flat(~infire);
        t_mask = tign_ground > tmax;
        tign_ground(t_mask) = tmax;
        tign_temp = imgaussfilt(tign_ground,1/6);
        %tign_temp = imgaussfilt(tign_ground,smooth_ground );
        %tign_temp = smooth_up(ps.red.fxlong,ps.red.fxlat,tign_ground,a,b);
        tign_ground(~infire) = tign_temp(~infire);
%         a = 1/10;%1/2-1/(2*i);
%         tign_ground(infire) = a*tign_ground(infire)+(1-a)*tign_temp(infire);
%         figure(73),scatter3(ps.red.fxlong(~infire),ps.red.fxlat(~infire),tign_ground(~infire))
%         pause(1/2)
        figure(159)
        contourf(ps.red.fxlong,ps.red.fxlat,tign_ground,20,'k'),hold on
        scatter(pts(1:5:end,2),pts(1:5:end,1),'*r'),hold off
        t_str = sprintf('Perimeter Shrinking \n Iteration %d',i);
        save_str = sprintf('perim_shrink_%d',i);
        xlabel('Lon'),ylabel('Lat'),title(t_str)
%         xl=[-121.5094 -121.2496];xlim(xl)
%         yl = [46.0367   46.2035];ylim(yl)
%         savefig(save_str);
%         saveas(gcf,[save_str '.png']);
        %figure(160),mesh(ps.red.fxlong,ps.red.fxlat,tign_ground);
        %pause(.5)
        %t_min(i) = min(tign_ground(:));
    end
end
tign_new = tign_ground;
%% end of post smoothing

tign_new = smooth_up(ps.red.fxlong,ps.red.fxlat,tign_new,a,b);



% %fix the ripple on the top
% r_mask = tign_new >=max(tign_new(:))-0.05;
% tign_new(r_mask) = max(tign_new(:));


figure(fig_num),mesh(ps.red.fxlong,ps.red.fxlat,tign_new-t0)
title(title_str)
hold on,scatter3(ps.grid_pts(:,2),ps.grid_pts(:,1),ps.points(:,3)-t0,'*r'),hold off

%plot the forecast for comparison
figure,mesh(ps.red.fxlong,ps.red.fxlat,ps.red.tign-t0)
title('Forecast')
hold on,scatter3(ps.grid_pts(:,2),ps.grid_pts(:,1),ps.points(:,3)-t0,'*r')

filled_contours(ps,tign_new)

end


