format compact 
A=eye(3);
% lexicographic unit cube
%    x y z
X0 = [0 0 0  %1
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
X = T*X0+rand(3,1)*ones(1,8);
X=X0;

u = [1, 2, 3]';
[K, F, Jg] = hexa(A,X,u);
eig(K)  % one eigenvalue zero

% test F
% for linear field mu with values V at nodes V*F should be -integral (grad mu) * u
gradmu=[4,5,6]; c=1;
V = gradmu*X + c;
vol=hexa_volume(X)
exact = - gradmu*u*vol; 
err_F = V*F - exact

% invariariant to random isometry
S=vrrotvec2mat(vrrotvec(randn(3,1),randn(3,1)))*diag(sign(randn(3,1))); 
[K0, ~, ~] = hexa(A,X,u);
[K1, ~, ~] = hexa(A,S*X,u);
err_iso=norm(K0-K1,'fro')

% test stretch
X3 = diag([1 1 2])*X0; 
[K3, ~, ~] = hexa(A,X3,u);
