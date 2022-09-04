program main
! need to include petsc files as well as use modules
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscvec.h>
use petscsys
use petscdm
use petscdmda
use petscvec
! add modules here
use double
use variables
use auxillary_conditions
use read_write
!use fdm
implicit none
! petsc datatypes
PetscInt           :: loc
PetscScalar        :: xp, fun
! datatypes
real(dp)           :: runtime
character(len=128) :: fmt1, fmt2, fmt3, fmt4
integer            :: values(8)
integer            :: i, iter

call PetscInitialize(PETSC_NULL_CHARACTER, ierr) 
if(ierr /= 0) stop "PETSc not intialized"

call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr); CHKERRQ(ierr)
call MPI_Comm_size(PETSC_COMM_WORLD, nproc, ierr); CHKERRQ(ierr)
runtime = MPI_Wtime()
fmt1 = "(1x, 'Date ', i2.2, '/', i2.2, '/', i4.4, '; time ', &
       & i2.2, ':', i2.2, ':', i2.2 )"

if(rank == 0)then
   call date_and_time(VALUES=values)
   write(*,fmt1) values(3), values(2), values(1), values(5:7)
   print*,'Running git version: ', VERSION
   print*,'Number of processes: ', nproc
endif

! Creates an object that will manage the communication of one-dimensional regular array data
! that is distributed across some processors. 
call DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, Nx, nvar, stencil_width, &
                  PETSC_NULL_INTEGER, da, ierr) ; CHKERRQ(ierr)
call DMSetFromOptions(da, ierr) ; CHKERRQ(ierr)
call DMSetUp(da, ierr); CHKERRQ(ierr)
call DMCreateGlobalVector(da, xg, ierr); CHKERRQ(ierr)
call DMCreateGlobalVector(da, ug, ierr); CHKERRQ(ierr)
! Returns the global (x,y,z) indices of the lower left corner and size of the local region, excluding ghost points.
call DMDAGetCorners(da, ibeg, PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, nloc, &
                    PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)
CHKERRQ(ierr)            
call VecDuplicate(ug, xg, ierr); CHKERRQ(ierr)
dx = (xmax - xmin) / dble(Nx)
do i = ibeg, ibeg+nloc-1
   xp = i*dx ; loc = i
   call initial_condition(xp, fun)
   call VecSetValues(xg, one, loc, xp, INSERT_VALUES, ierr); CHKERRQ(ierr)
   call VecSetValues(ug, one, loc, fun, INSERT_VALUES, ierr);   CHKERRQ(ierr)
enddo
call VecAssemblyBegin(xg, ierr);  CHKERRQ(ierr)
call VecAssemblyEnd(xg, ierr);    CHKERRQ(ierr)
call VecAssemblyBegin(ug, ierr);  CHKERRQ(ierr)
call VecAssemblyEnd(ug, ierr);    CHKERRQ(ierr)
call DMGetLocalVector(da, xl, ierr); CHKERRQ(ierr)
call DMGetLocalVector(da, ul, ierr); CHKERRQ(ierr)
call DMGlobalToLocalBegin(da, xg, INSERT_VALUES, xl, ierr); CHKERRQ(ierr)
call DMGlobalToLocalEnd(da, xg, INSERT_VALUES, xl, ierr); CHKERRQ(ierr)
call DMGlobalToLocalBegin(da, ug, INSERT_VALUES, ul, ierr); CHKERRQ(ierr)
call DMGlobalToLocalEnd(da, ug, INSERT_VALUES, ul, ierr); CHKERRQ(ierr)

time = 0.d0
dt = 0.1d0 ! cfl * dx / (speed_x + 1.d-13)
iter = 0

call save_solution(iter)

!do while (time .lt. final_time)
!  if (time + dt .gt. final_time) then
!  endif
!  call solve_lin_adv_diff()
!  time = time + dt
!  iter = iter + 1
!  print*, 'time = ', time
!enddo

! Print elapsed time
runtime = MPI_Wtime() - runtime
fmt4 = "(1x,'Iter  time min =',e20.12)"
if(rank == 0)then
   write(*,fmt4) runtime/60.0
endif

call DMDestroy(da, ierr); CHKERRQ(ierr)
call VecDestroy(xg, ierr); CHKERRQ(ierr)
call VecDestroy(ug, ierr); CHKERRQ(ierr)
call PetscFinalize(ierr); CHKERRA(ierr)

end program main
