module variables
   implicit none

   integer,parameter :: dp = kind(1.0d0)
   
   ! nproc = number of parallel processes
   ! rank  = index of this process (starts at 0)
   integer :: nproc, rank
   
   integer :: nx
   real    :: xmin, xmax, dx,  dt, dtp
   integer :: itmax
   integer :: itsave
   real    :: cfl
   real    :: final_time



   integer :: max_pit
   real    :: RESTOL, gamma=1.4
   integer :: schemetype
   integer :: WENO0 = 0, WENO1 = 1, WENO2 = 2, WENO3 = 3
   integer :: WENO5 = 5 , WENO4= 4
  

   integer :: nrk
   real    :: ark(3)
   real    :: gh, gl
   real :: M_PI = 4.0*atan(1.0), eps = 1.0e-12

   integer :: fileid_sol, fileid_omg

   integer :: fluxtype
   integer :: iroe=1, irusanov=2



   integer :: no=0, yes=1
   
   integer :: nvar=8
   integer :: var_rho = 1
   integer :: var_m1  = 2
   integer :: var_m2  = 3
   integer :: var_m3  = 4
   integer :: var_e   = 5
   integer :: var_b1  = 6
   integer :: var_b2  = 7
   integer :: var_b3  = 8

    ! WENO Constant
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   real    :: p107520=1.0/107520.0, p80640=1.0/80640.0,p8=1.0/8.0
   real    :: p228096=1.0/228096.0, p2112=1.0/2112.0
   real    :: p37440=1.0/37440.0, p720=1.0/720.0
   real    :: p5040=1.0/5040.0, p10=1.0/10.0, p126=1.0/126.0
   real    :: p1716=1.0/1716.0, p3=1.0/3.0
   real    :: p2002=1.0/2002.0, p2275=1.0/2275.0,p27=1.0/27.0
   real    :: p60918=1.0/60918.0, p93816426=1.0/93816426.0
   real    :: p1377684=1.0/1377684.0, p559521548066=1.0/559521548066.0
   real    :: p27582029244=1.0/27582029244.0, p443141066068272=1.0/443141066068272.0
   real    :: p560=1.0/560.0, p144=1.0/144.0
   real    :: p20=1.0/20.0, p12=1.0/12.0, p4=1.0/4.0, p40=1.0/40.0, p6=1.0/6.0
   real    :: p455=1.0/455.0, p49203=1.0/49203.0
   real    :: p90=1.0/90.0, p36=1.0/36.0
   real    :: p446=1.0/446.0, p336=1.0/336.0, p4480=1.0/4480.0, p24=1.0/24.0
   real    :: p120=1.0/120.0, p9=1.0/9.0

   real    :: p10080=1.0/10080.0, p216=1.0/216.0, p1584=1.0/1584.0
   real    :: p240=1.0/240.0, p56=1.0/56.0, p420=1.0/420.0, p60=1.0/60.0

 

end module variables
