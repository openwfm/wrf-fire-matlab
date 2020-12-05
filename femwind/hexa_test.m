format compact 
A=eye(3);
% lexicographic unit cube
%    x y z
X = [0 0 0
     1 0 0
     0 1 0
     1 1 0
     0 0 1
     1 0 1
     0 1 1
     1 1 1]';
 
u = [1, 2, 3]';
[K, F] = hexa(A,X,u)
eig(K)  % one eigenvalue zero

shp = alphaShape(X');
h = plot(shp);
vol=volume(shp)
% for linear field mu with values V at nodes V*F should be -integral (grad mu) * u
gradmu=[4,5,6]; c=1;
mu = @(x) gradmu*x + c;
V = mu(X);
ex = - gradmu*u*vol
err_F = V*F - ex
