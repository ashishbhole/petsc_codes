module read_write
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
use petscsys
use petscvec
use double
use variables
implicit none
contains

subroutine save_solution(iter)
implicit none
integer, intent(in)  :: iter
VecScatter     :: ctx
Vec            :: uall, xall
character(30) :: filename
PetscScalar, pointer :: xa(:), ua(:)
integer       :: i

!call VecView(ug, PETSC_VIEWER_STDOUT_WORLD, ierr); CHKERRQ(ierr)

call VecScatterCreateToZero(xg, ctx, xall, ierr); CHKERRQ(ierr)
call VecScatterCreateToZero(ug, ctx, uall, ierr); CHKERRQ(ierr)
call VecScatterBegin(ctx, ug, uall, INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
call VecScatterEnd(ctx, ug, uall, INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
call VecScatterBegin(ctx, xg, xall, INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
call VecScatterEnd(ctx, xg, xall, INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
call VecScatterDestroy(ctx, ierr); CHKERRQ(ierr)

!call VecView(uall, PETSC_VIEWER_STDOUT_WORLD, ierr); CHKERRQ(ierr)
if(rank==0)then
  call VecGetArrayF90(xall, xa, ierr); CHKERRQ(ierr)        
  call VecGetArrayF90(uall, ua, ierr); CHKERRQ(ierr)
  write(filename, '(a,i7.7,a)') 'solution_', iter, '.dat'
  open(10,file=trim(filename))
  write(10,*) '# xg, ug'
  do i = 1, Nx
     write(10, *) xa(i), ua(i)
  enddo  
  close(10)
  call VecRestoreArrayF90(xall, xa, ierr); CHKERRQ(ierr)  
  call VecRestoreArrayF90(uall, ua, ierr); CHKERRQ(ierr)  
endif
call VecDestroy(xall,ierr); CHKERRQ(ierr)
call VecDestroy(uall,ierr); CHKERRQ(ierr)

end subroutine save_solution

end module read_write
