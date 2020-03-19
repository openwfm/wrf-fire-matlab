function plot_conditions(n,v,C,D,desc)

figure, h = axes;
plot(C*v,'Color','k'),
hold on
nc1 = (n(1)-1)*n(2)*n(3);
line([nc1 nc1],get(h,'YLim'),'Color',[1 0 0],'LineWidth',2)
nc2 = nc1 + n(1)*(n(2)-1)*n(3);
line([nc2 nc2],get(h,'YLim'),'Color',[0 1 0],'LineWidth',2)
nc3 = nc2 + n(1)*n(2)*(n(3)-1);
line([nc3 nc3],get(h,'YLim'),'Color',[0 0 1],'LineWidth',2)
hold off
title(sprintf('Continuity conditions %s',desc))

figure,
plot(D*v)
title(sprintf('Divergence-free condition %s',desc))

end
