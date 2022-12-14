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

  ! constants
  real(dp)           :: pi = 4.d0*atan(1.d0), eps = 1.d-12
  PetscInt           :: one = 1, zero = 0

  PetscErrorCode     :: ierr
  PetscInt           :: rank, nproc
  
  DM                 :: da
  TS                 :: ts
  Vec                :: ug, ue
  TSType             :: time_scheme

  ! auxillary condition
  PetscReal          :: speed_x = 1.d0, speed_y = 1.d0
  real(dp)           :: amplitude = 1.d0, alpha = 64.d0, x_0 = 0.5d0, y_0= 0.5d0
  
  logical            :: petsc_ts
  character*64       :: space_disc
  integer            :: stencil_type         ! 1 for star, 2 for box
  PetscInt           :: stencil_width = 20

  PetscInt           :: Nx = 100, Ny = 100
  PetscReal          :: dx, dy, xmin=0.d0, xmax=1.d0, ymin=0.d0, ymax=1.d0
  PetscReal          :: dt, cfl, time=0.d0, final_time
  PetscInt           :: ibeg, jbeg, Nx_loc, Ny_loc
  PetscInt           :: ibeg_ghosted, Nx_loc_ghosted
  PetscInt           :: jbeg_ghosted, Ny_loc_ghosted
  integer            :: ist, ien, gist, gien ! gist : ghosted ist and so on
  integer            :: jst, jen, gjst, gjen

  real(dp)           :: err_l2
  integer            :: iter, itmax=100, itsave=10
  integer            :: igridformat, iascii = 1, ihdf5 = 2
  
  ! if required, construct a user defined data type here
  type tsdata
     real(dp)        :: dummy
  end type tsdata

end module variables
