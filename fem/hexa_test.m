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
% K = hexa(A,X)
hexa(A,X)