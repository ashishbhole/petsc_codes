module read_write
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdmda.h>
use petscsys
use petscvec
use petscdmda
use double
use variables
implicit none
contains

subroutine save_solution(iter, ug)
implicit none
PetscInt       :: iter
Vec            :: ug
DM             :: da
Vec            :: ul
Vec            :: uall
character(30) :: filename
PetscScalar, pointer :: ua(:,:)
integer       :: i, j
PetscViewer    viewer;

call TSGetDM(ts, da, ierr)
CHKERRQ(ierr)
call DMGetLocalVector(da, ul, ierr); CHKERRQ(ierr)
call DMGlobalToLocalBegin(da, ug, INSERT_VALUES, ul, ierr)
CHKERRQ(ierr)
call DMGlobalToLocalEnd(da, ug, INSERT_VALUES, ul, ierr)
CHKERRQ(ierr)
call DMDAVecGetArrayF90(da, ul, ua, ierr); CHKERRQ(ierr)
!call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr); CHKERRQ(ierr)
write(filename, '(a,i3.3,a,i3.3,a)') 'solution_', iter, '_', rank, '.plt'
open(10,file=trim(filename))
write(10,*) "TITLE = solution"
write(10,*) "VARIABLES = x, y, sol"
write(10,*) "ZONE STRANDID=1, SOLUTIONTIME=", time ,", I=", ien-ist+1, &
            ", J=", jen-jst+1, ", DATAPACKING=POINT"
do j = jst, jen
do i = ist, ien
   write(10, *) xmin+(i-1)*dx, ymin+(j-1)*dy, ua(i,j)
enddo
enddo
close(10)
call DMDAVecRestoreArrayF90(da, ul, ua, ierr); CHKERRQ(ierr)
call DMRestoreLocalVector(da, ul, ierr); CHKERRQ(ierr)
end subroutine save_solution

subroutine save_solution_hdf5(iter, da, ug)
implicit none
PetscInt       :: iter
DM             :: da
Vec            :: ug
PetscViewer    viewer;
character(30) :: filename

write(filename, '(a)') 'solution.h5'
call PetscViewerHDF5Open(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, viewer, ierr)
CHKERRQ(ierr)
call VecView(ug, viewer, ierr); CHKERRQ(ierr)
call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr);
end subroutine save_solution_hdf5

end module read_write
