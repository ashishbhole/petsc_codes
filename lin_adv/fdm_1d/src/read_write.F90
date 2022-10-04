module read_write
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
use petscsys
use petscvec
use double
use variables
implicit none
contains

subroutine save_solution(iter, ctx)
implicit none
PetscInt       :: iter
type(tsdata)   :: ctx
VecScatter     :: vsc
Vec            :: uall, xall
character(30) :: filename
PetscScalar, pointer :: xa(:), ua(:)
integer       :: i

call VecScatterCreateToZero(ctx%g%xg, vsc, xall, ierr)
CHKERRQ(ierr)
call VecScatterCreateToZero(ug, vsc, uall, ierr)
CHKERRQ(ierr)
call VecScatterBegin(vsc, ug, uall, INSERT_VALUES, SCATTER_FORWARD, ierr)
CHKERRQ(ierr)
call VecScatterEnd(vsc, ug, uall, INSERT_VALUES, SCATTER_FORWARD, ierr)
CHKERRQ(ierr)
call VecScatterBegin(vsc, ctx%g%xg, xall, INSERT_VALUES, SCATTER_FORWARD, ierr)
CHKERRQ(ierr)
call VecScatterEnd(vsc, ctx%g%xg, xall, INSERT_VALUES, SCATTER_FORWARD, ierr)
CHKERRQ(ierr)
call VecScatterDestroy(vsc, ierr); CHKERRQ(ierr)

!call VecView(uall, PETSC_VIEWER_STDOUT_WORLD, ierr); CHKERRQ(ierr)
if(rank==0)then
  call VecGetArrayF90(xall, xa, ierr); CHKERRQ(ierr)        
  call VecGetArrayF90(uall, ua, ierr); CHKERRQ(ierr)
  write(filename, '(a,i7.7,a)') 'solution_', iter, '.dat'
  open(10,file=trim(filename))
  write(10,*) '# xg, ug'
  do i = 1, ctx%g%Np
     write(10, *) xa(i), ua(i)
  enddo  
  close(10)
  call VecRestoreArrayF90(xall, xa, ierr); CHKERRQ(ierr)  
  call VecRestoreArrayF90(uall, ua, ierr); CHKERRQ(ierr)  
endif
call VecDestroy(xall,ierr); CHKERRQ(ierr)
call VecDestroy(uall,ierr); CHKERRQ(ierr)

end subroutine save_solution

subroutine set_and_braodcast_parameters(ctx)
implicit none
type(tsdata)   :: ctx
 
if(rank == 0)then

  select case(trim(space_disc))
  case('CD2')
    stencil_width = 3
  case('CD4')
    stencil_width = 5
  case('CD6')
    stencil_width = 7
  case('CD8')
    stencil_width = 9
  case('LELE')
    stencil_width = 20
  case default
    print*, 'please select a finite difference method in space'
    print*, 'Available options: CD2'
    stop
  end select

endif

call MPI_Barrier(PETSC_COMM_WORLD, ierr)

call MPI_Bcast(stencil_width, 1, MPI_int, 0, PETSC_COMM_WORLD, ierr)

end subroutine set_and_braodcast_parameters

subroutine log_parameters(ctx)
implicit none
type(tsdata)   :: ctx

! Check that floating points have correct precision
if (precision(pi) .ne. precision(pp)) then
  print*, 'Mismatch in float precision'
  stop
endif

print*,'*************About PETSc*************************'
      write(*,25) PETSC_VERSION_MAJOR, PETSC_VERSION_MINOR, &
                  PETSC_VERSION_SUBMINOR ! , PETSC_VERSION_PATCH
25    format(' PETSc Release Version: ',i1,'.',i2,'.',i1) ! ,' Patch:',i2)
print*,'****************************************************'

print*,'*************Discretization*************************'
print*, 'No of points = ', ctx%g%Np
print*, 'grid size = ', ctx%g%dx
print*, 'CFL number = ', ctx%g%cfl
print*, 'time step = ', ctx%g%dt
print*,'****************************************************'

print*,'*************Initial condition**********************'
print*, 'Amplitude of Gaussian Function= ', amplitude
print*, 'exponent = ', alpha
print*, 'centered around = ', x_0
print*,'****************************************************'

print*,'*************Finite difference methods**************'
print*, 'Use PETSc TS = ', petsc_ts 
print*, 'FDM in space = ', space_disc
if(.not. petsc_ts) print*, 'FDM in time = ', time_scheme
print*,'****************************************************'

print*,'*************Parallel computations******************'
print*, 'number of processes = ', nproc
print*, 'stencil width = ', stencil_width
print*,'****************************************************'

end subroutine log_parameters

end module read_write
