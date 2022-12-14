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
   PetscInt           :: stencil_width
   PetscReal          :: speed = 1.d0

   DM                 :: da
   TS                 :: ts
   TSType             :: time_scheme
   Vec                :: ug, ue
   
   ! auxillary condition
   real(dp) :: amplitude = 1.d0, alpha = 128.d0, x_0 = 0.4d0
   
   character*64  :: space_disc = 'CD6', time_disc ! 0 : explicit, 1 : implicit
   integer  :: nrk
   real(dp) :: ark(3)

   ! petsc functionalities
   logical  :: petsc_ts = .true.

   integer  :: ist, ien, gist, gien

   type grid
      PetscInt    :: Np = 100
      PetscReal   :: dx, xmin=0.d0, xmax=1.d0
      integer     :: iter, itmax=10000, itsave=10
      PetscReal   :: dt, cfl = 0.1, time, final_time=0.2d0      
      PetscInt    :: ibeg, nloc, ibeg_ghosted, nloc_ghosted   
  end type grid

   type tsdata
      type(grid)  :: g
   end type tsdata
   
   real(dp) :: err_l2

end module variables
