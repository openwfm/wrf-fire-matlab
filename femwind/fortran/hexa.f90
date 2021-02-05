subroutine hexa(A,X,u0,Kloc,Floc,Jg,iflags)
! purpose: create local stiffness matrix etc for hexa 3d element
! in:
!   A   coefficient matrix size 3x3, symmetric positive definite
!   X   nodes coordinates size 3x8, one each column is one node
!   u0   column vector of input wind at the center of the element
!   iflags  iflag(1)>0 compute Kloc, iflag(2)>0 compute Floc, iflag(3)>0 compute Jg  
! out:
!   Kloc   local stiffness matrix
!   Floc   local divergence load vector
!   Jg     gradient at center of function with values V is V'*Jg          

implicit none

!*** arguments
 
REAL, INTENT(IN):: A(3,3), X(3,8), u0(3)    ! fortran is not case sensitive
INTEGER, INTENT(IN)::iflag(3)
REAL, INTENT(OUT):: Kloc(8,8), Floc(8), Jg(3)

!*** local variables

!*** executable

! replace this by fortran code

!% basis functions on reference element [-1,1]^3
!Nb = 8;  % number of basis functions
!ib =[  % coordinates of basis functions
!    -1    -1    -1
!    -1    -1     1
!    -1     1    -1
!    -1     1     1
!     1    -1    -1
!     1    -1     1
!     1     1    -1
!     1     1     1];
!% the value of basis function k at x is
!% bf= @(k,x) (1+ib(k,1)*x(1))*(1+ib(k,2)*x(2))*(1+ib(k,3)*x(3))/8;
!%check_symmetry(A,'A',eps)
!% gaussian quadrature nodes
!g=0.5773502691896257;
!s = g*ib;
!Ng=Nb;  %  number of Gauss points
!s(Ng+1,:)=0; % extra point at center
!
!Kloc = zeros(Nb);
!Floc = zeros(Nb,1);
!for j=1:Ng+1
!    gradf = gradbfs(s(j,:));
!    Jx = X*gradf; % Jacobian at s
!    [q,r]=qr(Jx);
!    Jg   = (gradf/r)*q'; % gradf*inv(Jx)
!    adetJx = abs(prod(diag(r))); % det(Jx)
!    if j<=Ng % contribution to stiffness
!        K_at_s = Jg * A * Jg' * adetJx;
!        Kloc = Kloc + K_at_s;
!    else   % contribution to divergence load
!        % vol = adetJx * 8;
!        % TO DO: this is exact for a linearly deformed mesh but squashed is not.
!        % the decomposition in tetras used in hexa_volume will break
!        % non-planar faces. Compare the derminant and the average height method
!        % instead.
!        vol = hexa_volume(X);
!        Floc = Floc - Jg * u0 * vol;
!    end
!end
!%check_symmetry(Kloc,'Kloc',eps)
!end
    
end
