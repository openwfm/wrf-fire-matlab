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
[K, F] = hexa(A,10*X,u)
eig(K)