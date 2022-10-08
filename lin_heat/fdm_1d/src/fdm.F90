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
PetscScalar  :: ua(gist:gien), res(ist:ien)
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
case default
  print*, 'please select a finite difference method in space'
  print*, 'Available options: CD2'
  stop
end select

end subroutine finite_diffence_method


subroutine CD2(ua, res, ctx)
implicit none
type(tsdata) :: ctx
PetscScalar  :: ua(gist:gien), res(ist:ien)
! Local variables
integer      :: i
real(dp)     :: idx

idx = 1.d0/ctx%g%dx
do i = ist, ien
   res(i) = nu * (ua(i+1)-2.d0*ua(i)+ua(i-1))
enddo
res = res * idx * idx
end subroutine CD2

subroutine CD4(ua, res, ctx)
implicit none
type(tsdata) :: ctx
PetscScalar  :: ua(gist:gien), res(ist:ien)
! Local variables
integer        :: i
real(dp) :: idx, a_m2, a_m1, a_0, a_p1, a_p2

a_m2 = - 1.d0/12.d0
a_m1 =  4.d0/3.d0
a_0  = - 5.d0/2.d0
a_p1 =  a_m1
a_p2 =  a_m2
idx  = 1.d0/ctx%g%dx
do i = ist, ien
   res(i) = nu * ( a_m2*ua(i-2) + a_m1*ua(i-1) + a_0*ua(i) + &
                     a_p1*ua(i+1) + a_p2*ua(i+2) )
enddo
res = res * idx * idx
end subroutine CD4

subroutine CD6(ua, res, ctx)
implicit none
type(tsdata) :: ctx
PetscScalar  :: ua(gist:gien), res(ist:ien)
! Local variables
integer        :: i
real(dp) :: idx, a_m3, a_m2, a_m1, a_0, a_p1, a_p2, a_p3

a_m3 =   1.d0/90.d0
a_m2 = - 3.d0/20.d0
a_m1 =   3.d0/2.d0
a_0  = - 49.d0/18.d0
a_p1 =  a_m1
a_p2 =  a_m2
a_p3 =  a_m3
idx  = 1.d0/ctx%g%dx
do i = ist, ien
   res(i) = nu * ( a_m3*ua(i-3) + a_m2*ua(i-2) + a_m1*ua(i-1) + a_0*ua(i) + &
                   a_p1*ua(i+1) + a_p2*ua(i+2) + a_p3*ua(i+3) )
enddo
res = res * idx * idx
end subroutine CD6

subroutine CD8(ua, res, ctx)
implicit none
type(tsdata) :: ctx
PetscScalar  :: ua(gist:gien), res(ist:ien)
! Local variables
integer        :: i
real(dp) :: idx, a_m4, a_m3, a_m2, a_m1, a_0, a_p1, a_p2, a_p3, a_p4

a_m4 = - 1.d0/560.d0
a_m3 =   8.d0/315.d0
a_m2 = - 1.d0/5.d0
a_m1 =   8.d0/5.d0
a_0  = - 205.d0/72.d0
a_p1 =  a_m1
a_p2 =  a_m2
a_p3 =  a_m3
a_p4 =  a_m4
idx  = 1.d0/ctx%g%dx
do i = ist, ien
   res(i) = nu * ( a_m4*ua(i-4) + a_m3*ua(i-3) + a_m2*ua(i-2) + a_m1*ua(i-1) + a_0*ua(i) + &
                   a_p1*ua(i+1) + a_p2*ua(i+2) + a_p3*ua(i+3) + a_p4*ua(i+4) )
enddo
res = res * idx * idx
end subroutine CD8

end module fdm
