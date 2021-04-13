program restriction_test

use module_coarsening   
use module_io_matlab ! to read and write matrices as text files from matlab

implicit none

integer::                                         &
    ifds, ifde, kfds,kfde, jfds,jfde,             & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms,jfme,             &
    ifts, ifte, kfts,kfte, jfts,jfte

integer::                             &
    ifcds, ifcde, kfcds,kfcde, jfcds,jfcde,       & ! coarse grid domain
    ifcms, ifcme, kfcms,kfcme, jfcms,jfcme,       & ! coarse grid dimensions
    ifcts, ifcte, kfcts,kfcte, jfcts,jfcte          ! coarse grid tile

real, pointer, dimension(:,:,:):: u_m, uc_m, X_m, Y_m, Z_m, cl_z_m, hcz_m ! to read from files
real, pointer, dimension(:,:,:):: u, uc, X, Y, Z  ! to pass on
integer, pointer:: cl_z(:)  
integer:: n(3),nc(3),cr_x,cr_y, nwrap = 0

! read matrices, X Y Z already in the ikj ordering
call read_array(X_m,'X')
call read_array(Y_m,'Y')
call read_array(Z_m,'Z')
call read_array(u_m,'u')
call read_array(uc_m,'uc')  ! for dimensions only
call read_array(cl_z_m,'cl_z')
call read_array(hcz_m,'hcz')

n = shape(X_m)     ! fine mesh size
nc = shape(uc_m)   ! coarse mesh size

call set_indices(nc,nwrap,                  &    ! coarse mesh bounds
    ifcds, ifcde, kfcds,kfcde, jfcds,jfcde, & 
    ifcms, ifcme, kfcms,kfcme, jfcms,jfcme, &
    ifcts, ifcte, kfcts,kfcte, jfcts,jfcte)


! set coarsening indices
if(nc(2).ne.size(cl_z_m))call crash('coarse index array size')
allocate(cl_z(kfcts:kfcte))
cl_z = cl_z_m(:,1,1);
cr_x = hcz_m(1,1,1)
cr_y = hcz_m(2,1,1)

call set_indices(n,nwrap,                         &     ! coarse mesh bounds
    ifds, ifde, kfds,kfde, jfds,jfde,             & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms,jfme,             &
    ifts, ifte, kfts,kfte, jfts,jfte)

! allocate fine input and copy to tile 
allocate(X(ifms:ifme,kfms:kfme,jfms:jfme))
X(ifts:ifte,kfts:kfte,jfts:jfte)=X_m
allocate(Y(ifms:ifme,kfms:kfme,jfms:jfme))
Y(ifts:ifte,kfts:kfte,jfts:jfte)=Y_m
allocate(Z(ifms:ifme,kfms:kfme,jfms:jfme))
Z(ifts:ifte,kfts:kfte,jfts:jfte)=Z_m
allocate(u(ifms:ifme,kfms:kfme,jfms:jfme))
u(ifts:ifte,kfts:kfte,jfts:jfte)=u_m

! allocate output
allocate(uc(ifcms:ifcme,kfcms:kfcme,jfcms:jfcme))

write(*,*)'calling restriction'
call restriction(   &
    ifds, ifde, kfds,kfde, jfds, jfde,            & ! fire grid dimensions
    ifms, ifme, kfms,kfme, jfms, jfme,            &
    ifts, ifte, kfts, kfte, jfts,jfte,            &
    ifcds, ifcde, kfcds,kfcde, jfcds,jfcde,       & ! coarse grid domain
    ifcms, ifcme, kfcms,kfcme, jfcms,jfcme,       & ! coarse grid dimensions
    ifcts, ifcte, kfcts,kfcte, jfcts,jfcte,       & ! coarse grid tile
    uc,u,cr_x,cr_y,cl_z,X,Y,Z)

! allocate coarse output and copy
allocate(uc_m(nc(1),nc(2),nc(3)))
uc_m=uc(ifcts:ifcte,kfcts:kfcte,jfcts:jfcte)
call write_array(uc_m,'uc')

end program restriction_test
