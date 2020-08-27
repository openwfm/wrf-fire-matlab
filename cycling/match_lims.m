function match_lims(a,b)
%function copies xlim,ylim from figure(a) to figure(b)

figure(a)
xl = xlim;
yl = ylim;

figure(b)
xlim(xl)
ylim(yl)

end
