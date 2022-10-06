module auxillary_conditions
#include <petsc/finclude/petscsys.h>
use petscsys
use double
use variables
implicit none
contains

subroutine initial_condition(xp, fun)
implicit none
integer :: i
PetscScalar  :: xp, fun
fun = amplitude * exp(-alpha*(xp - x_0)**2)
end subroutine initial_condition

subroutine exact_solution(xp, fun, ctx)
implicit none
integer :: i
PetscScalar  :: xp, fun
type(tsdata) :: ctx
PetscReal    :: k1, delta

k1 = 4.0d0 * nu * ctx%g%final_time
delta = (1.d0/k1 + alpha)
fun = ( 1.d0/ dsqrt(4.0d0 * nu * ctx%g%final_time * delta) ) &
          * exp( - ((xp-x_0)**2 / k1) * (1.0d0-1.0d0/(delta * k1)) )
end subroutine exact_solution

end module auxillary_conditions
