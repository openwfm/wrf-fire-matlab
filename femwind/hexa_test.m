format compact 
A=eye(3);
% lexicographic unit cube
%    x y z
X = [0 0 0  %1
     1 0 0  %2
     0 1 0  %3
     1 1 0  %4
     0 0 1  %5
     1 0 1  %6
     0 1 1  %7
     1 1 1  %8
     ]';

%   7-----8
%  /|    /|
% 5-----6 |
% | |   | |
% | 3---|-4
% |/    |/
% 1-----2 

% linear transformation
T = rand(3);
% T = magic(3);
X = T*X;

u = [1, 2, 3]';
[K, F] = hexa(A,X,u)
eig(K)  % one eigenvalue zero

vol=hexa_volume(X)
% for linear field mu with values V at nodes V*F should be -integral (grad mu) * u
gradmu=[4,5,6]; c=1;
mu = @(x) gradmu*x + c;
V = mu(X);
ex = - gradmu*u*vol
err_F = V*F - ex
