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
!case('CD6')
!  call CD6(ua, res, ctx)
!case('CD8')  
!  call CD8(ua, res, ctx)
case('LELE')
  call LELE(ua, res, ctx)
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

idx = 0.5d0/ctx%g%dx
do i = ist, ien
   res(i) = - speed * (ua(i+1)-ua(i-1))
enddo
res = res * idx
end subroutine CD2

subroutine CD4(ua, res, ctx)
implicit none
type(tsdata) :: ctx
PetscScalar  :: ua(gist:gien), res(ist:ien)
! Local variables
integer        :: i
real(dp) :: idx, a_m2, a_m1, a_p1, a_p2

a_m2 =  1.d0/12.d0
a_m1 = -2.d0/3.d0
a_p1 = - a_m1
a_p2 = - a_m2
idx  = 1.d0/ctx%g%dx
do i = ist, ien
   res(i) = - speed * ( a_m2*ua(i-2) + a_m1*ua(i-1) + &
                        a_p1*ua(i+1) + a_p2*ua(i+2) )
enddo
res = res * idx
end subroutine CD4

subroutine LELE(ua, res, ctx)
implicit none
type(tsdata) :: ctx
PetscScalar  :: ua(gist:gien), res(ist:ien)
! Local variables
integer      :: i
real(dp)     :: udia(gist:gien), dia(gist:gien)
real(dp)     :: ldia(gist:gien), rhs(gist:gien)
real(dp)     :: res_tmp(gist:gien)
real(dp) :: a_m1, a_0, a_p1, idx
real(dp) :: b_m2, b_m1, b_0, b_p1, b_p2

udia = 0.0 ; dia = 0.0 ; ldia = 0.0; rhs = 0.0; res_tmp = 0.0

!interior nodes
a_m1 = 1.0/3.0
a_0  = 1.0
a_p1 = a_m1
b_m2 = -1.0/36.0
b_m1 = -14.0/18.0
b_0  = 0.0
b_p1 = - b_m1
b_p2 = - b_m2

idx = 1.0/ctx%g%dx
ldia(gist) = 0.0 ; dia(gist) = 1.0 ; udia(gist) = 0.0 
rhs(gist) = - speed * (ua(gist+1) - ua(gist)) * idx
ldia(gist+1) = 0.0 ; dia(gist+1)  = 1.0 ; udia(gist+1) = 0.0
rhs(gist+1) = - speed * (0.5*ua(gist+2) - 0.5* ua(gist)) * idx
do i = gist+2, gien-2
   ldia(i) = a_m1; dia(i) = a_0; udia(i) = a_p1
   rhs(i) = - speed * ( b_m2 * ua(i-2) + b_m1 * ua(i-1) + b_0 * ua(i) + &
                        b_p1 * ua(i+1) + b_p2 * ua(i+2) ) * idx
enddo
ldia(gien-1) = 0.0 ; dia(gien-1) = 1.0 ; udia(gien-1) = 0.0
rhs(gien-1) = - speed * (0.5*ua(gien) - 0.5*ua(gien-2)) * idx
ldia(gien) = 0.0 ; dia(gien)  = 1.0 ; udia(gien) = 0.0
rhs(gien)  = - speed * (ua(gien) - ua(gien-1)) * idx

call tdma(ldia(gist+1:gien), dia(gist:gien), udia(gist:gien-1), & 
           rhs(gist:gien), res_tmp(gist:gien), gien-gist+1)

res(ist:ien) = res_tmp(ist:ien)

end subroutine LELE

subroutine tdma(a3, b3, c3, d3, xa, N)
implicit none
integer, intent(in):: N
real(dp), intent(in), dimension(2:N)  :: a3
real(dp), intent(in), dimension(1:N-1):: c3
real(dp), intent(in), dimension(1:N)  :: b3, d3
PetscScalar, intent(out), dimension(1:N) :: xa
real(dp), allocatable, dimension(:) :: beta1, gamma1
integer :: i

allocate(beta1(1:N),gamma1(1:N))

beta1(1)=b3(1)
gamma1(1)=d3(1)/b3(1)
do i = 2, N
  beta1(i)  =   b3(i) - a3(i) * c3(i-1)/beta1(i-1)
  gamma1(i) = ( d3(i) - a3(i)*gamma1(i-1) ) / beta1(i)
enddo

xa(N) = gamma1(N)
do i = N-1, 1, -1
  xa(i) = gamma1(i) - c3(i) * xa(i+1)/beta1(i)
enddo

deallocate(beta1, gamma1)

end subroutine TDMA

end module fdm
