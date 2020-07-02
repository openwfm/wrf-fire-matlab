function tign_new = squish(ps)

tign_new = ps.red.tign;
tign = ps.red.tign;
idx = ps.idx;
mask = ones(size(tign));
%change back to z = tign
% if exist('ps.tn','var')
%     figure(68),contour3(ps.red.fxlong,ps.red.fxlat,ps.tn-ps.red.start_datenum,20)
% else
%     figure(68),contour3(ps.red.fxlong,ps.red.fxlat,tign-ps.red.start_datenum,20)
% end
% hold on
%lons = ps.red.fxlong;
%lats = ps.red.fxlat;


%combine l2 time with detection points fixed to grid
pts = ps.grid_pts;
pts(:,3) = ps.points(:,3);


%just squish the ignition point
for k = 1:12
for i = 1:length(ps.paths)
    p = ps.paths(i).p;
    plot3(pts(p,2),pts(p,1),pts(p,3)-ps.red.start_datenum,'r')
    %plot3(pts(p,2),pts(p,1),tign(idx(p,1),idx(p,2))-ps.red.start_datenum,'g')
    for j = 1:length(p)
        p_i = idx(p(j),1)+round(randn);
        p_j = idx(p(j),2)+round(randn);
        tign_new(p_i-round(rand):p_i+round(rand),p_j-round(rand):p_j+round(rand)) = pts(p(j),3)-0.25*rand;
%        mask(idx(p(j),1),idx(p(j),2)) = mask(idx(p(j),1),idx(p(j),2))-1/length(ps.paths);
%         holes = length(tign(:))-sum(mask(:));
%         fprintf('%d holes in mask \n',holes);
        if j > 1
            new_i = uint8(round((idx(p(j),1)+idx(p(j-1),1))/2));
            new_j = uint8(round((idx(p(j),2)+idx(p(j-1),2))/2));
            new_t =  0.5*(pts(p(j),3)+pts(p(j-1),3));
            tign_new(new_i-round(rand):new_i+round(rand),new_j-round(rand):new_j+round(rand)) = new_t-0.25*rand;
            %mask(new_i,new_j)=0;
        end
        for k = 1:1
       % tign_new = rlx_shp(tign_new,1/2,mask);
        end
    end
end
tign_new = rlx_shp(tign_new,1/2,mask);
end
%tign_new = [];
figure,mesh(ps.red.fxlong,ps.red.fxlat,tign_new)
end

function new_mask = expand_mask(mask,border)
    [n,m] = size(mask);
    [ix,jx] = find(mask(border:n-border,border:m-border)==0);
    for i = 1:length(ix)
        mask(ix(i)-border:ix(i)+border,jx(i)-border:jx(i)+border)=1-1/border;
    end  
    new_mask = mask;
end