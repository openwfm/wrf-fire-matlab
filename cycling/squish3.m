function tign_new = squish3(ps)
%% ps is struct with paths, red, graph,distances, etc.
%% ps = graph_dets(w)


%combine l2 time with detection points fixed to grid
pts = ps.grid_pts;
pts(:,3) = ps.points(:,3);
%forecast
tign = ps.red.tign;
t0 = min(tign(:));

tign2 = squish2(ps);
tign2 = smooth_up(ps.red.fxlong,ps.red.fxlat,tign2);
tign = max(tign,tign2);
idx = ps.idx;
fig_num = 23;
pts_length = length(ps.grid_pts);
%max time to look at detection data
max_l2_time = max(max(pts(:,3)));


tign_new = ps.red.tign;
%tign_flat is constant fire arrival time
tign_flat=ps.red.end_datenum*ones(size(tign));
title_str = 'Analysis';

%make matrix of forecast tign time for the points in the graph
t_times = zeros(pts_length,1);
for i = 1:pts_length
    t_times(i) = ps.red.tign(ps.idx(i,1),ps.idx(i,2));
end

%set up data likelihood spline
[p_like_spline,~,~] = make_spline(100,1000);
%test
% t = 4*(-10:20)*3600;
% ts = p_like_spline(t);
% figure,plot(t,exp(ts))
%make vector of data likelihoods

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
rm = 0;

% random multiplier, keep the same
% perturbs points downward in time to
rt = 0.2;
% weight for tign_new

%alhpa blends  estimate of tign at a point with old estimate
% new_estimate = alph*current_setimate + (1-alpha)*old_estimate
% for data assimilation, alpha will be computed from exp(likelihood)
alpha = 0.5;
%constant for smooth in rlx_shp
alpha_2 = 0.5; %smaller alph_2 ==> smoother
%number of loops to run
smoothings = 10;
for k = 1:smoothings
    figure(fig_num),mesh(ps.red.fxlong,ps.red.fxlat,tign_new)
    title(title_str)
    hold on,scatter3(ps.grid_pts(:,2),ps.grid_pts(:,1),ps.points(:,3),'*r'),hold off
    %pause(3/k)
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
            alpha = alpha_vect(p(j));
            tign_new(p_i-round(rm*rand):p_i+round(rm*rand),p_j-round(rm*rand):p_j+round(rm*rand)) = alpha*tign_new(p_i,p_j) + (1-alpha)*pts(p(j),3)-rt*rand;
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
    patch = 4;
    
    %smooth the tign
    tign_new(tign_new < t0) = t0;
    tign_new = smooth_up(ps.red.fxlong,ps.red.fxlat,tign_new);
    tign_new = rlx_shp(tign_new,alpha_2,patch);
%     tign_flat = rlx_shp(tign_flat,alpha_2,patch);
    
    %collect information about tign at the detection points
    for i = 1:pts_length
        t_times(i) = tign_new(ps.idx(i,1),ps.idx(i,2));
        flat_times(i) = tign_flat(ps.idx(i,1),ps.idx(i,2));
    end
    norms(k,1) = norm(t_times-ps.points(:,3),2);
    % blend flat start with forecast start for analysis
    
    %tign_flat = rlx_shp(tign_flat,alpha_2,patch);
    norm_flat = norm(flat_times-ps.points(:,3),2);
    norm_tign_new = norms(k,1);
    if (norm_flat > norm_tign_new)
        %fprintf('blending flat start with forecast \n')
        %tign_new = 0.5*(tign_new+tign_flat);
    end
    
    %only do norm for times before final detection time
    time_mask = tign_new < pts(end,3);  %max(max(pts(:,3)));
    norms(k,2) = norm(tign_new(time_mask)-tign_old(time_mask));
    %     figure(fig_num+3);plot(norms(1:k,1));
    %     tstr = sprintf('Norm of difference between successive \n TIGN after each interpolation');
    %     title(tstr)
    %     figure(fig_num+4);plot(1:k,norms(1:k,2))
    %     tstr = sprintf('Norm of difference between times of \n detections and TIGN at detectionon locations');
    %     title(tstr);
    %     xlabel('Iterations of Interpolation')
    temp_var(k) = min(tign_new(:))-ps.red.start_datenum;
    figure(fig_num+5);plot(1:k,temp_var(1:k)*24),title('Change in ignition time'),xlabel('iteration'),ylabel('hours')
    fprintf('Loop %d complete norm of diff = %f \n', k,norms(k))
    if k > 2 && norms(k,1) > norm(k-1,1) && norms(k,1) > norms(k-2,1)
        fprintf('graph norm increase \n')
        %break
    end
end
%tign_new = [];

figure,mesh(ps.red.fxlong,ps.red.fxlat,ps.red.tign)
title('Forecast')
hold on,scatter3(ps.grid_pts(:,2),ps.grid_pts(:,1),ps.points(:,3),'*r')

end

function new_mask = expand_mask(mask,border)
[n,m] = size(mask);
[ix,jx] = find(mask(border:n-border,border:m-border)==0);
for i = 1:length(ix)
    mask(ix(i)-border:ix(i)+border,jx(i)-border:jx(i)+border)=1-1/border;
end
new_mask = mask;
end


