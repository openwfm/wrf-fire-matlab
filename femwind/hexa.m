function [Kloc,Floc,Jg]=hexa(A,X,u0)
% create local stiffness matrix for hexa 3d element
% in:
%   A   coefficient matrix size 3x3, symmetric positive definite
%   X   nodes coordinates size 3x8, one each column is one node
%   u0   column vector of input wind at the center of the element
% out:
%   Kloc   local stiffness matrix
%   Floc   local divergence load vector
%   Jg     gradient at center of function with values V is V'*Jg          

% basis functions on reference element [-1,1]^3
Nb = 8;  % number of basis functions
ib =[  % coordinates of basis functions
    -1    -1    -1
     1    -1    -1
    -1     1    -1
     1     1    -1
    -1    -1     1
     1    -1     1
    -1     1     1
     1     1     1];
ibt=ib';
% the value of basis function k at x is
% bf= @(k,x) (1+ib(k,1)*x(1))*(1+ib(k,2)*x(2))*(1+ib(k,3)*x(3))/8;
%check_symmetry(A,'A',eps)
% gaussian quadrature nodes
g=0.5773502691896257;
s = g*ib;
Ng=Nb;  %  number of Gauss points
s(Ng+1,:)=0; % extra point at center

Kloc = zeros(Nb);
Floc = zeros(Nb,1);
for j=1:Ng+1
    gradf = gradbfs(s(j,:));
    Jx = X*gradf; % Jacobian at s
    [q,r]=qr(Jx);
    Jg   = (gradf/r)*q'; % gradf*inv(Jx)
    adetJx = abs(prod(diag(r))); % det(Jx)
    if j<=Ng % contribution to stiffness
        K_at_s = Jg * A * Jg' * adetJx;
        Kloc = Kloc + K_at_s;
    else   % contribution to divergence load
        % find volume 
        % rectangle base in x1 and x2 times average height in x3
        % we really should be passing dx(1) dx(2) as arguments not compute them
        dx(1) = X(1,:)*ib(:,1)/4; % average length of sides in x1 direction
        dx(2) = X(2,:)*ib(:,2)/4; % average length of sides in x2 direction
        dx(3) = X(3,:)*ib(:,3)/4; % average length of sides in x3 direction
        vola=dx(1)*dx(2)*dx(3);
        % determinant at midpoint times volume of reference element
        % exact if element is linear image of reference element
        vold = adetJx*8;   
        % volume by decompositon into tetras
        % exact if each side has corners in a plane
        % but sides have kinks if not
        volt = hexa_volume(X);   
        vol_err=(abs([vola,vold]-volt))/volt
        vol=vold;
        Floc = Floc - Jg * u0 * vol;
    end
end
%check_symmetry(Kloc,'Kloc',eps)
end
    
