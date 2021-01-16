n=10;
la=linspace(0,3,n);
ld=linspace(0,3,n);
a=10.^la;
d=10.^ld;
for i=1:n
    for j=1:n
        C=hexa_coupling([1,1,a(i)],[1,1,d(j)]);
        cv(i,j)=C(1,5);
        ch(i,j)=C(1,3);
        laa(i,j)=la(i);
        ldd(i,j)=ld(j);
        aa(i,j)=a(i);
        dd(i,j)=d(j);
    end
end
figure(1)
mesh(laa,ldd,cv)
hold on
contour(laa,ldd,cv,[-1:0.1:1])
hold off
grid on
title('Vertical coupling')
xlabel('log vertical coefficient')
ylabel('log vertical size')

figure(2)
mesh(laa,ldd,ch)
hold on
contour(laa,ldd,ch,[-1:0.1:1])
hold off
grid on
title('Horizontal coupling')
xlabel('log vertical coefficient')
ylabel('log vertical size')