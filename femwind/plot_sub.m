function plot_sub(X,K,i)
% plot_sub(X,K,i)
% plot submatrix of K for element i=[i1 i2 i3] from the mesh X{1} X{2} X{3}
p=hexa_sub(X,K,i);
hexa_plot(p.X,p.X,p.S2);

