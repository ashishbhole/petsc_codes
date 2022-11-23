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

! Compute rhs function in yg' = r(t,yg), where yg and r are global vectors
subroutine RHSFunction(ts, time, yg, r, ctx, ierr)
  implicit none
  TS             :: ts
  PetscReal      :: time
  Vec            :: yg
  Vec            :: r
  type(tsdata)   :: ctx
  DM             :: da
  Vec            :: localU   
  PetscErrorCode :: ierr
  PetscScalar, pointer :: u(:,:), res(:,:)
  
  call TSGetDM(ts, da, ierr); CHKERRQ(ierr)
  call DMGetLocalVector(da, localU, ierr); CHKERRQ(ierr)
  call DMGlobalToLocalBegin(da, yg, INSERT_VALUES, localU, ierr); CHKERRQ(ierr)
  call DMGlobalToLocalEnd(da, yg, INSERT_VALUES, localU, ierr); CHKERRQ(ierr)
  call DMDAVecGetArrayReadF90(da, localU, u, ierr); CHKERRQ(ierr)
  call DMDAVecGetArrayF90(da, r, res, ierr); CHKERRQ(ierr)
  call finite_diffence_method(u(gist:gien,gjst:gjen), res(ist:ien,jst:jen))
  call DMDAVecRestoreArrayReadF90(da, localU, u, ierr); CHKERRQ(ierr)
  call DMDAVecRestoreArrayF90(da, r, res, ierr); CHKERRQ(ierr)
  call DMRestoreLocalVector(da, localU ,ierr); CHKERRQ(ierr)
end subroutine RHSFunction

! can be used to compute time step, save solution and monitor estimates
subroutine Monitor(ts, step, time, u, ctx, ierr)
  implicit none
  TS             :: ts
  PetscInt       :: step
  PetscReal      :: time, u_l2
  Vec            :: u
  type(tsdata)   :: ctx
  PetscErrorCode :: ierr
  
  call VecNorm(u, NORM_2, u_l2, ierr); CHKERRQ(ierr)
  if(rank == 0)then
    write(*,'(a, i4, 3F8.4)') 'iter, dt, time, l2_norm = ', step, dt, time, u_l2
  endif
  if( mod(step, itsave) == 0 )then
    if(iascii == igridformat)then
      call save_solution(step, ug)
    elseif(ihdf5 == igridformat)then
      call save_solution_hdf5(step, ug)
    endif
  endif
end subroutine Monitor

end module solve
