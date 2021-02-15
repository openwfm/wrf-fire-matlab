function [mx,my,s] = moe_2d(w1,w2)
%w1 is wrfout from forecast
%w2 is wrfout with "ground truth" fire arrival time - optional if perimeter
%              is used
perim_use = 0;
if ~exist('w2','var')
    [fire_name,save_name,prefix,perim] = fire_choice()
    perim_use = 1;
end
[m,n] = size(w1.fxlong);
%perim_use = input_num('Compare to Perimeter? Yes = 1',1,1)

if perim_use == 1

    
    % load data if it has been processes already for another comparison
    load_perim = input_num('Load perimeter data? Yes = 1',1,0);
    if load_perim == 1 && exist('perim_moe.mat','file')
        load perim_moe.mat
    else
        perim_points = input_num('Number of perim points to use?',200)
        perim_struct = perim2gran(5000,perim);
        for i = 1:length(perim_struct)
            fprintf('%d : %s \n',i,perim_struct(i).file)
        end
        perim_choice = input_num('Which perimeter to use ?',3);
        %interpolate the perimeter to the grid
        perim_points = fixpoints2grid(w1,[perim_struct(perim_choice).lat',perim_struct(perim_choice).lon']);
        perim_points = unique(perim_points,'rows');
        %make polygon comparison
       
        [x,y] = meshgrid(1:n,1:m);
        in_perim = inpolygon(y,x,perim_points(:,1),perim_points(:,2));
        in_perim = inpolygon(y,x,x(in_perim),y(in_perim));
        fprintf('Choose default bounds 2\n')
        r = subset_domain(w1,0);
        tign_g = r.tign_g; %datenum2time(w1.tign_g,r);
        perim_time = perim_struct(perim_choice).time;
        t_g = datenum2time(perim_time,r);
        t_g = r.end_time;
        
        save perim_moe.mat perim_points in_perim x y m n perim_struct perim_choice r perim_time t_g tign_g
    end
    

    perim_flat = zeros(m,n);
    perim_flat(in_perim) = 1;
    perim_area = sum(perim_flat(:));
    
else
    t_g = min(max(w1.tign_g(:)),max(w2.tign_g(:)))-6*3600;
    tign_g = w1.tign_g;
    %area of ground truth
    in_perim = w2.tign_g<=t_g;
    
    perim_flat = zeros(m,n);
    perim_flat(in_perim') = 1;
    perim_area = sum(perim_flat(:));
    
    
end %if perim_use...

%start comparison with w.tign_g

w1_mask = tign_g < t_g;
w1_flat = zeros(m,n);;
w1_flat(w1_mask') = 1;
w1_area = sum(w1_flat(:));

%find overlap
overlap = perim_flat+w1_flat;
overlap(overlap<1.5) = 0;
overlap(overlap>0) = 1;
overlap = logical(overlap);

false_neg = perim_flat-w1_flat;
false_neg(false_neg < 0) = 0;
false_neg = logical(false_neg);
false_pos = w1_flat - perim_flat;
false_pos(false_pos < 0) = 0;
false_pos = logical(false_pos);

%MOE
mx = 1 - sum(false_neg(:))/perim_area;
my = 1 - sum(false_pos(:))/w1_area;
%s --- sorenson index
s = 2*sum(overlap(:))/(w1_area+perim_area);

plot_on = 1;
if plot_on == 1
    figure(129)
    close(gcf)
    figure(129),scatter(w1.fxlong(false_neg'),w1.fxlat(false_neg'),'r*')
    hold on,scatter(w1.fxlong(false_pos'),w1.fxlat(false_pos'),'k*'),
    scatter(w1.fxlong(overlap'),w1.fxlat(overlap'),'g*')
    legend('False Negative','False Positive','Overlap')
    tstr = sprintf('MOE = (%f,%f) \n Sorenson = %f',mx,my,s);
    title(tstr);
end

end