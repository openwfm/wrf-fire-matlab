function [K,F]=hexa(A,X,u0)
% create local stiffness matrix for hexa 3d element
% in:
%   A   coefficient matrix size 3x3, symmetric positive definite
%   X   nodes coordinates size 3x8, one each column is one node
%   u0   column vector of input wind at the center of the element
% out:
%   K   local stiffness matrix

% basis functions on reference element [-1,1]^3
Nb = 8;  % number of basis functions
k=0; ib=zeros(Nb,3);
for i1=-1:2:1,for i2=-1:2:1,for i3=-1:2:1
    k=k+1; 
    ib(k,:)=[i1,i2,i3];  %  coordinates of node k in the reference element 
end,end,end
% the value of basis function k at x
bf= @(k,x) (1+ib(k,1)*x(1))*(1+ib(k,2)*x(2))*(1+ib(k,3)*x(3))/8;
% bd{i}(k,x) is the derivative wrt x(i) of basis function k at x
gradbf = @(k,x) ...  % gradient of basis function k at x
     [ib(k,1)*(1+ib(k,2)*x(2))*(1+ib(k,3)*x(3)),... % d/d(x(1))
      (1+ib(k,1)*x(1))*ib(k,2)*(1+ib(k,3)*x(3)),... % d/d(x(2))
      (1+ib(k,1)*x(1))*(1+ib(k,2)*x(2))*ib(k,3)]/8; % d/d(x(3))

% gaussian quadrature nodes
g=0.5773502691896257;
s = g*ib;
Ng=Nb;  %  number of Gauss points
s(Ng+1,:)=0; % extra point at center

% gradient on Gauss points - precompute
for j=1:Ng+1
    for k=1:Nb
        gradfs(k,:,j)= gradbf(k,s(j,:));
    end
end

% changes with X
K = zeros(Nb);
F = zeros(Nb,1);
for j=1:Ng+1
    gradf = gradfs(:,:,j);
    Jx = X*gradf; % Jacobian at s
    [q,r]=qr(Jx);
    Jg   = (gradf/r)*q'; % gradf/Jx
    adetJx = abs(prod(diag(r))); % det(Jx)
    if j<=Ng % contribution to stiffness
        K_at_s = Jg * A * Jg' * adetJx;
        K = K + K_at_s;
    else   % contribiution to divergence load
        F = F + Jg * u0 * adetJx * 8;
    end
end
    
