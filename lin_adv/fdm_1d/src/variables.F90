module variables
   use double     
   implicit none

   real(dp) :: pi = 4.d0*atan(1.d0), eps = 1.d-12, tolerance
   integer  :: nproc, rank ! number of processes and their index
   
   integer  :: nx
   real(dp) :: xmin, xmax, dx, dt, cfl, final_time
   integer  :: itmax, itsave

   integer  :: space_disc, time_disc ! 0 : explicit, 1 : implicit
   integer  :: nrk
   real(dp) :: ark(3)

   integer :: nvar=1
   integer :: var_u = 1

end module variables
