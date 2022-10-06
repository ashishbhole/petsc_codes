module fdm
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
use petscsys
use petscvec
use double
use variables
implicit none
contains

subroutine finite_diffence_method(ua, res)
implicit none
PetscScalar :: ua(gist:gien,gjst:gjen), res(ist:ien,jst:jen)
! Local variables

res = 0.d0

select case(trim(space_disc))
case('CD2')
  call CD2(ua, res)
!case('LELE')
!  call LELE(ua, res)
case default
  print*, 'please select a finite difference method in space'
  print*, 'Available options: CD2'
  stop
end select

end subroutine finite_diffence_method


subroutine CD2(ua, res)
implicit none
PetscScalar :: ua(gist:gien,gjst:gjen), res(ist:ien,jst:jen)
! Local variables
integer        :: i, j
real(dp) :: idx, idy

idx = 0.5d0/dx; idy = 0.5d0/dy
do j = jst, jen
do i = ist, ien
   res(i,j) = - speed_x * (ua(i+1,j)-ua(i-1,j)) * idx &
            + - speed_y * (ua(i,j+1)-ua(i,j-1)) * idy
enddo
enddo

end subroutine CD2

!subroutine LELE(ua, res)
!implicit none
!PetscScalar,pointer  :: res(:,:), ua(:,:)
! Local variables
!integer        :: i, ist, ien
!real(dp), allocatable, dimension(:) :: udia, dia, ldia, rhs
!real(dp) :: a_m1, a_0, a_p1, idx
!real(dp) :: b_m2, b_m1, b_0, b_p1, b_p2


!ist = ctx%g%ibeg-stencil_width
!ien = ctx%g%ibeg+ctx%g%nloc+stencil_width
!
!allocate(udia(ist:ien), dia(ist:ien), ldia(ist:ien), rhs(ist:ien))
!udia = 0.0 ; dia = 0.0 ; ldia = 0.0; rhs = 0.0
!
!!interior nodes
!a_m1 = 1.0/3.0
!a_0  = 1.0
!a_p1 = a_m1
!b_m2 = -1.0/36.0
!b_m1 = -14.0/18.0
!b_0  = 0.0
!b_p1 = - b_m1
!b_p2 = - b_m2
!
!idx = 1.0/ctx%g%dx
!ldia(ist) = 0.0 ; dia(ist) = 1.0 ; udia(ist) = 0.0 
!rhs(ist) = - speed * (ua(ist+1) - ua(ist)) * idx
!ldia(ist+1) = 0.0 ; dia(ist+1)  = 1.0 ; udia(ist+1) = 0.0
!rhs(ist+1) = - speed * (0.5*ua(ist+2) - 0.5* ua(ist)) * idx
!do i = ist+2, ien-2
!   ldia(i) = a_m1; dia(i) = a_0; udia(i) = a_p1
!   rhs(i) = - speed * ( b_m2 * ua(i-2) + b_m1 * ua(i-1) + b_0 * ua(i) + &
!                        b_p1 * ua(i+1) + b_p2 * ua(i+2) ) * idx
!enddo
!ldia(ien-1) = 0.0 ; dia(ien-1) = 1.0 ; udia(ien-1) = 0.0
!rhs(ien-1) = - speed * (0.5*ua(ien) - 0.5*ua(ien-2)) * idx
!ldia(ien) = 0.0 ; dia(ien)  = 1.0 ; udia(ien) = 0.0
!rhs(ien)  = - speed * (ua(ien) - ua(ien-1)) * idx
!
!call tdma(ldia(ist+1:ien), dia(ist:ien), udia(ist:ien-1), rhs(ist:ien), res(ist:ien), ien-ist+1)
!
!deallocate(udia, dia, ldia)
!
!end subroutine LELE

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
