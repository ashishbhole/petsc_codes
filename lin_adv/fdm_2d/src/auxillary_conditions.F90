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

end module auxillary_conditions