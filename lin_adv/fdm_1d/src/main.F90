program main
! need to include petsc files as well as use modules
#include <petsc/finclude/petscsys.h>
use petscsys
! add modules here
use double
use variables
!use fvm
implicit none
! petsc datatypes
PetscErrorCode     :: ierr
PetscReal, allocatable :: cons, prim, res

! datatypes
integer :: i
real(dp)           :: runtime
character(len=128) :: fmt1, fmt2, fmt3, fmt4
integer            :: values(8)

call PetscInitialize(PETSC_NULL_CHARACTER, ierr) 
if(ierr /= 0) stop "PETSc not intialized"

call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
CHKERRA(ierr)
call MPI_Comm_size(PETSC_COMM_WORLD, nproc, ierr)
CHKERRA(ierr)
runtime = MPI_Wtime()
fmt1 = "(1x, 'Date ', i2.2, '/', i2.2, '/', i4.4, '; time ', &
       & i2.2, ':', i2.2, ':', i2.2 )"

if(rank == 0)then
   call date_and_time(VALUES=values)
   write(*,fmt1) values(3), values(2), values(1), values(5:7)
   print*,'Running git version: ', VERSION
   print*,'Number of processes: ', nproc
endif

! Print elapsed time
runtime = MPI_Wtime() - runtime
fmt4 = "(1x,'Iter  time min =',e20.12)"
if(rank == 0)then
   write(*,fmt4) runtime/60.0
endif

call PetscFinalize(ierr); CHKERRA(ierr)

end program main
