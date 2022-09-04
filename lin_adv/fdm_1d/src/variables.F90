module variables
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscvec.h>
use petscsys
use petscdm
use petscdmda
use petscvec

   use double
   implicit none

   real(dp) :: pi = 4.d0*atan(1.d0), eps = 1.d-12, tolerance

   PetscErrorCode     :: ierr
   PetscInt           :: rank, nproc, nvar = 1, Nx = 100, one = 1, zero = 0
   PetscInt           :: stencil_width = 3
   PetscReal          :: dx, xmin=0.d0, xmax=1.d0, dt, cfl, time, final_time=1.d0, speed_x

   DM                 :: da
   PetscInt           :: ibeg, nloc
   Vec                :: xg, xl, ug, ul
   
   integer  :: itmax, itsave
   integer  :: space_disc, time_disc ! 0 : explicit, 1 : implicit
   integer  :: nrk
   real(dp) :: ark(3)

   !integer :: nvar=1
   integer :: var_u = 1

   ! auxillary condition
   real(dp) :: amplitude = 1.d0, alpha = 64.d0, x_0 = 0.5d0

end module variables
