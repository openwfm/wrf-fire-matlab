t=[1 1 1];
% t=sym('t',[1,3])
s =[
    -1    -1    -1
    -1    -1     1
    -1     1    -1
    -1     1     1
     1    -1    -1
     1    -1     1
     1     1    -1
     1     1     1];
s = s .* (ones(8,1)*t)
x=sym('x',[3,1]);
b=sym('b',[8,1]);
g=sym('g',[8,3]);
for k=1:8
    b(k) = (1+x(1)/s(k,1))*(1+x(2)/s(k,2))*(1+x(3)/s(k,3));
    for l=1:3
        g(k,l)=diff(b(k),x(l));
    end
end
K=int(int(int(g*g.',x(1),-t(1),t(1)),x(2),-t(2),t(2)),x(3),-t(3),t(3))