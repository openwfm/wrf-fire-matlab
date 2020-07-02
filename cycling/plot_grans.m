function plot_grans(g)

n = length(g);
b = g(1).time;
figure
for i = 2:n
    t1 = g(i-1).time - b;
    t2 = g(i).time -b;
    if mod(i,2)
        c = 'green';
    else
        c = 'red';
    end
    x = [0 n n 0];
    y = [t1 t1 t2 t2];
    hold on
    patch(x,y,c)

end