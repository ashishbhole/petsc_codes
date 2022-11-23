module auxillary_conditions
#include <petsc/finclude/petscsys.h>
use petscsys
use double
use variables
implicit none
contains

subroutine initial_condition(xp, yp, fun)
  implicit none
  integer :: i
  PetscScalar  :: xp, yp, fun
  fun = amplitude * exp(-alpha*((xp - x_0)**2 + (yp - y_0)**2))
end subroutine initial_condition

subroutine exact_solution(xp, yp, fun)
  implicit none
  integer :: i
  PetscScalar  :: xp, yp, fun
  fun = amplitude * exp(-alpha*(  ((xp-speed_x*final_time) - x_0)**2 &
                                + ((yp-speed_y*final_time) - y_0)**2 ) )
end subroutine exact_solution

end module auxillary_conditions
