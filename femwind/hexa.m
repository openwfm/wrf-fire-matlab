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
ib =[  % coordinates of basis functions
    -1    -1    -1
    -1    -1     1
    -1     1    -1
    -1     1     1
     1    -1    -1
     1    -1     1
     1     1    -1
     1     1     1];
% the value of basis function k at x is
% bf= @(k,x) (1+ib(k,1)*x(1))*(1+ib(k,2)*x(2))*(1+ib(k,3)*x(3))/8;

% gaussian quadrature nodes
g=0.5773502691896257;
s = g*ib;
Ng=Nb;  %  number of Gauss points
s(Ng+1,:)=0; % extra point at center

K = zeros(Nb);
F = zeros(Nb,1);
for j=1:Ng+1
    gradf = gradbfs(s(j,:));
    Jx = X*gradf; % Jacobian at s
    [q,r]=qr(Jx);
    Jg   = (gradf/r)*q'; % gradf*inv(Jx)
    adetJx = abs(prod(diag(r))); % det(Jx)
    if j<=Ng % contribution to stiffness
        K_at_s = Jg * A * Jg' * adetJx;
        K = K + K_at_s;
    else   % contribution to divergence load
        F = F - Jg * u0 * adetJx * 8;
    end
end
end
    
