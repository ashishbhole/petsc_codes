module solve
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscts.h>
#include <petsc/finclude/petscdmda.h>
use petscsys
use petscts
use petscdmda
use double
use variables
use fdm
use read_write
implicit none
contains

! Compute rhs function in y' = f(t,y)
subroutine RHSFunction(ts, time, ug, r, ctx, ierr)
#include <petsc/finclude/petscts.h>
#include <petsc/finclude/petscdmda.h>
use petscts
use petscdmda
implicit none
TS             :: ts
DM             :: da
PetscReal      :: time
Vec            :: ug
Vec            :: r
type(tsdata)   :: ctx
Vec            :: localU   
PetscErrorCode :: ierr
PetscReal      :: val
PetscInt       :: i, loc
PetscScalar, pointer :: u(:), res(:)

call TSGetDM(ts, da, ierr)
CHKERRQ(ierr)
call DMGetLocalVector(da, localU, ierr)
CHKERRQ(ierr)
call DMGlobalToLocalBegin(da, ug, INSERT_VALUES, localU, ierr)
CHKERRQ(ierr)
call DMGlobalToLocalEnd(da, ug, INSERT_VALUES, localU, ierr)
CHKERRQ(ierr)
call DMDAVecGetArrayReadF90(da, localU, u, ierr)
CHKERRQ(ierr)
call DMDAVecGetArrayF90(da, r, res, ierr)
CHKERRQ(ierr)
call DMDAGetCorners(da, ctx%g%ibeg, PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, ctx%g%nloc, &
                 PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)
CHKERRQ(ierr)
call finite_diffence_method(u, res, ctx)
call DMDAVecRestoreArrayReadF90(da, localU, u, ierr)
CHKERRQ(ierr)
call DMDAVecRestoreArrayF90(da, r, res, ierr)
CHKERRQ(ierr)
call DMRestoreLocalVector(da, localU ,ierr)
CHKERRQ(ierr)
end subroutine RHSFunction

! This subroutine is called after every time step
! Note: This function seems to be called before time stepping starts.
subroutine Monitor(ts, step, time, u, ctx, ierr)
use petscts
implicit none
TS             :: ts
PetscInt       :: step
PetscReal      :: time
Vec            :: u
type(tsdata)   :: ctx
PetscErrorCode :: ierr

if(rank == 0) print*, 'iter, dt, time = ', step, ctx%g%dt, time
if( mod(step, ctx%g%itsave) == 0 ) call save_solution(step, ctx)

end subroutine Monitor

end module solve
