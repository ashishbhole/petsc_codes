module variables
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscts.h>
use petscsys
use petscdm
use petscdmda
use petscvec
use petscts

   use double
   implicit none

   real(dp) :: pi = 4.d0*atan(1.d0), eps = 1.d-12, tolerance
   PetscReal      :: pp

   PetscErrorCode     :: ierr
   PetscInt           :: rank, nproc, one = 1, zero = 0
   PetscInt           :: stencil_width = 3
   PetscReal          :: speed

   DM                 :: da
   TS                 :: ts
   TSType             :: time_scheme
   Vec                :: ug
   
   ! auxillary condition
   real(dp) :: amplitude = 1.d0, alpha = 64.d0, x_0 = 0.5d0
   
   character*64  :: space_disc = 'CD2', time_disc ! 0 : explicit, 1 : implicit
   integer  :: nrk
   real(dp) :: ark(3)

   ! petsc functionalities
   logical  :: petsc_ts = .true.

   type grid
      PetscInt    :: Np = 100
      PetscReal   :: dx, xmin=0.d0, xmax=1.d0
      integer     :: iter, itmax=100, itsave=10
      PetscReal   :: dt, cfl, time, final_time=1.d0      
      Vec         :: xg
      PetscInt    :: ibeg, nloc   
  end type grid

   type tsdata
      type(grid)  :: g
   end type tsdata
   
end module variables
