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
character(30) :: filename
PetscScalar, pointer :: ua(:,:)
integer       :: i, j
PetscViewer    viewer;

call TSGetDM(ts, da, ierr); CHKERRQ(ierr)
call DMGetLocalVector(da, ul, ierr); CHKERRQ(ierr)
call DMGlobalToLocalBegin(da, ug, INSERT_VALUES, ul, ierr); CHKERRQ(ierr)
call DMGlobalToLocalEnd(da, ug, INSERT_VALUES, ul, ierr); CHKERRQ(ierr)
call DMDAVecGetArrayF90(da, ul, ua, ierr); CHKERRQ(ierr)
write(filename, '(a,i3.3,a,i3.3,a)') 'solution_', rank, '_', iter, '.plt'
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
call VecDestroy(ul,ierr); CHKERRQ(ierr)

end subroutine save_solution

subroutine save_solution_hdf5(iter, da, ug)
implicit none
PetscInt       :: iter
DM             :: da
Vec            :: ug
PetscViewer    viewer;
character(30) :: filename

write(filename, '(a)') 'solution.h5'
call PetscViewerHDF5Open(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, viewer, ierr); CHKERRQ(ierr)
call VecView(ug, viewer, ierr); CHKERRQ(ierr)
call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

end subroutine save_solution_hdf5

subroutine read_and_set_parameters()
implicit none
PetscBool :: set
character(len=64)  :: string
character(len=20)  :: fmt1, fmt2, fmt3

   if(rank==0) print*,'Reading parameters from param.in'

   cfl = 0.0d0 ! Default cfl number
   call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-cfl',cfl,set,ierr); CHKERRQ(ierr)

   dt = 0.0d0 ! Default time step
   call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-dt',dt,set,ierr); CHKERRQ(ierr)


   itmax = 10000000 ! Default: large value
   call PetscOptionsGetInt (PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-itmax',itmax,set,ierr); CHKERRQ(ierr)

   itsave = 10000000 ! Default: large value
   call PetscOptionsGetInt (PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-itsave',itsave,set,ierr); CHKERRQ(ierr)
    
   ! Format of grid files
   igridformat = iascii ! Default: is ascii
   string = ''
   call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                              '-grid_format',string,set,ierr); CHKERRQ(ierr)
   if(trim(string) == 'hdf5')then
      igridformat = ihdf5
   else
      igridformat = iascii     
   endif

   select case(trim(space_disc))
   case('CD2')
     stencil_width = 1; stencil_type = 1
   case('CD4')
     stencil_width = 2; stencil_type = 1
   case('CD6')
     stencil_width = 3; stencil_type = 1
   case('LELE')
     stencil_width = 10; stencil_type = 1
   case default
     print*, 'please select a finite difference method in space'
     print*, 'Available options: CD2, CD4, CD6'
     stop
   end select

   call MPI_Bcast(stencil_width, 1, MPI_int, 0, PETSC_COMM_WORLD, ierr)

   ! Check options only on root process.
   ! TODO: After this, abort is called on on root process, need to fix this.
   if(rank > 0) return

   ! Print paramaters to screen
   fmt1 = "(5x,a10,i14)"
   fmt2 = "(5x,a10,e14.4)"
   fmt3 = "(5x,a10,L)"
   write(*,fmt2) 'tfinal   =', final_time
   write(*,fmt2) 'cfl      =', cfl
   write(*,fmt2) 'dt       =', dt
   write(*,fmt1) 'itmax    =', itmax
   write(*,fmt1) 'itsave   =', itsave
   write(*,fmt3) 'tscheme  =', petsc_ts

   if(cfl < zero .and. dt < zero)then
      write(*,*)'Both cfl and dt cannot be mzero'
   endif

   if(cfl > zero .and. dt > zero)then
      write(*,*)'Both cfl and dt cannot be specified'
      write(*,*)'Give one and set other to mzero'
   endif
   
end subroutine read_and_set_parameters

end module read_write
