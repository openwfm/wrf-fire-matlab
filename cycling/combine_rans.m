function c = combine_rans(n,w)

c = ran_cone(w);
for i = 1:n
    c1 = ran_cone(w);
    a1 = sum(sum(c1<max(c1(:))));
    a = sum(sum(c<max(c(:))));
    c = a/(a+a1)*c+a1/(a+a1)*c1;
    %c = max(c,c1);
    %c = 1/2*(c+c1);
end

%c = imgaussfilt(c,1);
c = smooth_up(c);
% figure(123),mesh(c);
%figure(124),contour(c,20,'k')

t_min = min(w.tign_g(:));
t_max = max(w.tign_g(:));
c_min = min(c(:));
c_max = max(c(:));

% days = round((t_max-t_min)/24/3600);
% stepper = (c_max-c_min)/days;
% for i = 1:days
%    c_min = min(c(:));
%    c_max = max(c(:));
%    msk =  c>c_min+i*stepper;
%    figure(217),scatter(w.fxlong(msk),w.fxlat(msk))
%    pause(3)
%    sum(msk(:));
%    m_rand = 1+(-1)^i/2
%    %m_rand = 1/10*(1+sin((i+1)*pi*rand/2))^2;
%    c(~msk)= m_rand/10*c(~msk);
% end

c_min = min(c(:));
c_max = max(c(:));

c = t_min + (t_max-t_min)/(c_max-c_min)*(c-c_min);
%figure(123),mesh(w.fxlong,w.fxlat,c);
figure(124),contour(w.fxlong,w.fxlat,c,20,'k')

%quick_mesh(c)


end%function
