function C=hexa_coupling(a,d)
% hexa_coupling(a,d)
% example: hexa_coupling([1,1,100],[1,1,1])

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
X=diag(d)*X0;
A=diag(a);
u = [1, 2, 3]';
[K, F, Jg] = hexa(A,X,u);
D=diag(diag(K).^(-1/2));
C=D*K*D;
end