function K=hexa(A,X)
% create local stiffness matrix for hexa 3d element
% in:
%   A   coefficient matrix size 3x3, symmetric positive definite
%   X   nodes coordinates size 3x8, one each column is one node
% out:
%   K   local stiffness matrix

% basis functions on reference element [-1,1]^3
Nb = 8;  % number of basis functions
k=0; ib=zeros(Nb,3);
for i1=-1:2:1,for i2=-1:2:1,for i3=-1:2:1
    k=k+1; ib(k,:)=-[i1,i2,i3];
    % the value of basis function k at x is
    bf= @(k,x) (1+ib(k,1)*x(1))*(1+ib(k,2)*x(2))*(1+ib(k,3)*x(3))/8;
    % bd{i}(k,x) is the derivative wrt x(i) of basis function k at x
    bd{1}= @(k,x) ib(k,1)*(1+ib(k,2)*x(2))*(1+ib(k,3)*x(3))/8; % d/d(x(1))
    bd{2}= @(k,x) (1+ib(k,1)*x(1))*ib(k,2)*(1+ib(k,3)*x(3))/8; % d/d(x(2))
    bd{3}= @(k,x) (1+ib(k,1)*x(1))*(1+ib(k,2)*x(2))*ib(k,3)/8; % d/d(x(3))
end,end,end

% gaussian quadrature nodes
g=0.5773502691896257;
s = g*ib;
Ng=Nb;  %  number of Gauss points 

% gradient on Gauss points - precompute
for j=1:Ng
    for k=1:Nb
        gradfs(k,:,j)=[bd{1}(k,s(j,:)),bd{2}(k,s(j,:)),bd{3}(k,s(j,:))];
    end
end

% changes with X
K = zeros(Nb);
for j=1:Ng
    gradf = gradfs(:,:,j);
    Jx = X*gradf; % Jacobian at s
    Jg   = gradf/Jx;
    K_at_s = Jg * A * Jg' * abs(det(Jx));
    K = K + K_at_s;
end
    
