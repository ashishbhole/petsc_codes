module fdm
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
use petscsys
use petscvec
use double
use variables
implicit none
contains

subroutine finite_diffence_method(ua, res, ctx)
implicit none
type(tsdata)         :: ctx
PetscScalar,pointer  :: res(:), ua(:)
! Local variables
integer        :: i

res = 0.d0

select case(trim(space_disc))
case('CD2')
  call CD2(ua, res, ctx)
case('CD4')  
  call CD4(ua, res, ctx)
case('CD6')
  call CD6(ua, res, ctx)
case('CD8')  
  call CD8(ua, res, ctx)
case('LELE')
  call LELE(ua, res, ctx)  
case default
  print*, 'please select explicit finite difference method in space'
  print*, 'Available options: CD2'
  stop
end select

end subroutine finite_diffence_method


subroutine CD2(ua, res, ctx)
implicit none
type(tsdata)         :: ctx
PetscScalar,pointer  :: res(:), ua(:)
! Local variables
integer        :: i

res(1) = - 0.5d0 * speed * (ua(2)-ua(ctx%g%Np-1)) / ctx%g%dx
do i = ctx%g%ibeg, ctx%g%ibeg+ctx%g%nloc
   res(i) = - 0.5d0 * speed * (ua(i+1)-ua(i-1)) / ctx%g%dx
enddo
res(ctx%g%Np) = - 0.5d0 * speed * (ua(2)-ua(ctx%g%Np-1)) / ctx%g%dx

end subroutine CD2

subroutine CD4(ua, res, ctx)
implicit none
type(tsdata)         :: ctx
PetscScalar,pointer  :: res(:), ua(:)
! Local variables
integer        :: i

end subroutine CD4

subroutine CD6(ua, res, ctx)
implicit none
type(tsdata)         :: ctx
PetscScalar,pointer  :: res(:), ua(:)
! Local variables
integer        :: i

end subroutine CD6

subroutine CD8(ua, res, ctx)
implicit none
type(tsdata)         :: ctx
PetscScalar,pointer  :: res(:), ua(:)
! Local variables
integer        :: i

end subroutine CD8

subroutine LELE(ua, res, ctx)
implicit none
type(tsdata)         :: ctx
PetscScalar,pointer  :: res(:), ua(:)
! Local variables
integer        :: i

end subroutine LELE

end module fdm
