module module_hexa   ! testing only

contains

subroutine hexa(A,X,u0,Kloc,Floc,Jg,iflags)
! purpose: create local stiffness matrix etc for hexa 3d element
! in:
!   A   coefficient matrix size 3x3, symmetric positive definite
!   X   nodes coordinates size 3x8, one each column is one node
!   u0   column vector of input wind at the center of the element
!   iflags  iflags(1)>0 compute Kloc, iflags(2)>0 compute Floc, iflags(3)>0 compute Jg  
! out:
!   Kloc   local stiffness matrix
!   Floc   local divergence load vector
!   Jg     gradient at center of function with values V is V'*Jg          

implicit none

!*** arguments

real, intent(in):: A(3,3), X(3,8), u0(3)    ! fortran is not case sensitive
integer, intent(in)::iflags(3)
real, intent(out):: Kloc(8,8), Floc(8), Jg(3)

!*** local variables
real, parameter :: Nb = 8
real, parameter :: Ng
real, parameter :: tmp
real, parameter :: g = 0.5773502691896257
real, dimension(8, 3) :: ib
real, dimension(8, 3) :: s
real, dimension(3,3) :: Jx
real, dimension(3,3) :: qdash
real, dimension(3,3) :: Q
real, dimension(3,3) :: R
real, dimension(8,3) :: gradf
real, dimension(8, 3) :: Jg
real, dimension(8, 3) :: Jg_tran
real, dimension(8, 3) :: Jg_tmp
real, dimension(8,8) :: Kloc
real, dimension(8,8) :: K_at_s
real, dimension(3,3) :: Q_tran
real, dimension(3,3) :: R_inv

!*** executable

! temporary, assign someting 
! to prevent compiler warning about unassigned variables
Kloc = 0 
Floc = 0
Jg = 0
s = 0

ib = reshape((/-1,-1,-1,-1,1,1,1,1,-1,-1,1,1,-1,-1,1,1,-1,1,-1,1,-1,1,-1,1/),shape(ib))
Ng = Nb

!current state of fortran code. still neds a lot of work to find mor efficient coding scheme
!this will most likely break since it is incomplete

do i = 1,9
    do j = 1,3
        s(i,j) = g*ib(i,j)
    end do
end do

do i = 1,9  !need to turn these into subrountines

        !gradf piece
    do j = 1,8
        do k = 1,3
            gradf(j,k) = ((1+ib(j,1)*s(i,1))*(1+ib(j,2)*s(i,2))*(1+ib(j,3)*s(i,3)))/8
        end do
    end do

        !Jacobian at s
    do j = 1,3
        do m = 1,3
            tmp = 0
            do n = 1,8
               tmp = tmp + X(m,n)*gradf(n,m)
                end do
                Jx(j,k) = tmp
        end do
    end do

        !QR Decoposition (long way)        
        !Get first column of Q
    tmp = 0
    do j = 1,3
        qdash(j,1) = Jx(j,1)
        tmp = tmp + qdash(j,1)*qdash(j,1)
    end do

    R(1,1) = sqrt(tmp)

    do j = 1,3
        Q(j,1) = qdash(j,1)/R(1,1)
        R(1,2) = R(1,2) + Q(j,1)*Jx(j,2)
    end do

        !Get second column of Q
    tmp = 0
    do j =1,3
        qdash(j,2) = Jx(j,2)-R(1,2)*Q(j,1)
        tmp = tmp + qdash(j,2)*qdash(j,2)
    end do

    R(2,2) = sqrt(tmp)

    do j =1,3
        Q(j,2) = qdash(j,2)/R(2,2)
        R(1,3) = R(1,3) + Q(j,1)*Jx(j,3)
        R(2,3) = R(2,3) + Q(j,2)*Jx(j,3)
    end do

        !Get second column of Q
    tmp = 0
    do j =1,3
        qdash(j,3) = Jx(j,3) - R(1,3)*Q(j,1) - R(2,3)*Q(j,2)
        tmp = tmp + qdash(j,3)*qdash(j,3)
    end do

    R(3,3) = sqrt(tmp)

    do j = 1,3
        Q(j,3) = qdash(j,3)/R(3,3)
    end do
        
        !det of Jx
    detJx = abs(R(1,1)*R(2,2)*R(3,3))

        !Jg piece
    do j = 1,3
        do k = 1,3
            Q_tran(k,j) = Q(j,k)
            R_inv(k,j) = R(j,k)/(detJx)
        end do
    end do

    do j=1,8
        tmp = 0
        do k = 1,3
            do m = 1,3
                tmp = tmp+gradf(j,n)*R_inv(n,m)
            end do
        end do
        Jg_tmp(j,k) = tmp
    end do

    do j=1,8
        tmp = 0
        do k = 1,3
            do m = 1,3
                tmp = tmp+Jg_tmp(j,n)*Q_trans(n,m)
            end do
        end do
        Jg(j,k) = tmp
    end do

    if i <= Ng

        !insert K stuff here
        !K_at_s stuff here
        do j = 1,8
            do k = 1,3
                Jg_tran(k,j) = Jg(j,k)*detJx
            end do
        end do

        do j = 1,3
            do k = 1,8
                do m =1,3
                    tmp = tmp + A(j,m)*Jg_trans(m,k)
                end do
            end do
            A(j,k) = tmp
        end do

        do j = 1,8
            do m = 1,8
                tmp = 0
                do n = 1,3
                    tmp = tmp + Jg(m,n)*A(n,m)
                end do
                K_at_s(j,k) = tmp
            end do
        end do 


        !Kloc stuff here
        do j = 1,8
            do k=1,8
                Kloc(j,k) = Kloc(j,k)-K_at_s(j,k)
            end do
        end do

    else
       
        !fill in what used to be hexa_volume piece
        !Floc piece
        do j = 1,8
            tmp = 0
            do m = 1,3
                tmp = tmp+Jg(j,k)*u0(k,1)
            end do
            tmp_mat(j,1) = tmp*vol
        end do

        do j = 1,8
            Floc(j,1) = Floc(j,1) - tmp_mat(j,1)
        end do
    end if
end do


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
    
end subroutine hexa

end module module_hexa
