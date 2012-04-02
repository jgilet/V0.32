  
  ! This fil reads parameters from the namelist, and calculate the paramters 
  ! needed for model integration.
  
  !                                                                       
  !!     parameters:                                                       
  !     NI,NJ,NK are the no. of points in each of the 3 co-ordintate direc
  !     L, DL are the characteristic length and depth scales.             
  !     R is the radius of the earth = 6371 km.                           
  !     NEQ is the number of equations and the number or unknowns         
  !     EPS is the Rosby number                                           
  !     AH and AV are the eddy diffusivities in the horizontal and vertica
  !     lexp is a constant in the equations denoting the length scale     
  !     LEN= 10**(4+lexp),  0<=lexp<=1                                    
  !     dtime is the time step                                            
  
  !     np is the number of points which makes up the salinity,temp,dens v
  !     depth data curve                                                  
  
  !     the h eqn. kaphinv is its inverse.                                
  !     set kappah=1/2.  kaphinv=1/kappah     
  
  INTEGER :: step,nsteps,save_steps,n,i,j,k,nblock,ik,inew,iold
  integer :: frame_int,ngraph2d,m,budget_int,pickup_int 
  !      REAL(kind=rc_kind) :: fdiv,p(maxout),coslat(0:NJ+1),                  
  !      REAL(kind=rc_kind) :: fdiv,p(maxout),colati(0:NJ+1),
  REAL(kind=rc_kind) ::  fdiv,ctrdiv,hsum                                            
  REAL(kind=rc_kind) ::  tarray(2), tim 
  INTEGER :: lexp              
  REAL(kind=rc_kind) ::  pcorr(maxout)
  !real*8, allocatable :: pcorr(:)    
  character(len=10) :: stepchar
  !=========================================
  !      read configuration file
  NAMELIST /PARAM/ nsteps,save_steps,dtf,ngraph2d,frame_int,&
                   pickup_int,pickup_step,fnhhy,dx,dy,total_depth,Kx_TS,&
                   Ky_TS_bdy,Ky_TS_int,Kx_m,Ky_m,dirout,phi0deg
#ifdef particle
  NAMELIST /traj/ ini_particle_time,parti_file_num,NPR,pcx,pcy,pcz,pcr
#endif
  READ (*,PARAM) !from the namelist file
#ifdef particle
  READ (*,traj) !from the namelist file
#endif

  !calculate the basic parameter variables
  budget_int= 10 
  !     define the parameters                                             
  !=      frame_int= 100  read in gulf.in                                 
  !      EPS= 0.1d0 defined in header                                     
  delta= DL/LEN 
  delinv= LEN/DL
  !     set qpr to zero for hydrostatic, delta for nonhydrostatic         
  !      qpr= 0.d0                                                        
  qpr= fnhhy*delta 
  !      qprinv= 1.d0/qpr not to be used    
  kappah= 0.65d0 
  kaphinv= 1.d0/kappah 
  !=      kappah= 0.5d0                                                   
  !=      kaphinv= 2.d0                                                   
  lambda= LEN/AL 
  !     beta= 1.d0/(EPS*EPS*delta)  - orignal value                       
  beta= 1.d0/(EPS*EPS*delta) 
  !      beta=100.d0                                                      
  UL= FPAR *LEN*EPS 
  !     UL= F*LEN*EPS                                                     
  P1= R0*UL*UL*EPS**(-1) 
  !     P1= RU^2/EPS                                                      
  HL= P1/(R0*10.d0) 
  !     R0*G*HL= P1                                                       
  !     HDL= 1.d-4 or 1.d-3                                               
  HDL= HL/DL 
  TL= LEN/UL 
  WL= EPS*delta*UL 

  WRITE(*,"(A,E14.4)")  "# fnhhy   = ", fnhhy
  WRITE(*,"(A,E14.4)")  "# EPS     = ",EPS
  WRITE(*,"(A,E14.4)")  "# delta   = ",delta
  WRITE(*,"(A,E14.4)")  "# qpr     = ", qpr
  WRITE(*,"(A,E14.4)")  "# lambda  = ",lambda
  WRITE(*,"(A,E14.4)")  "# beta    = ", beta
  WRITE(*,"(A,E14.4)") "# dtf     = ", dtf
  WRITE(*,"(A,E14.4)")  "# WL      = ", WL
  WRITE(*,"(A,E14.4)")  "# R0      = ", R0
  WRITE(*,"(A,E14.4)")    "# UL      = ", UL
  WRITE(*,"(A,E14.4)")     "# P1      = ", P1
  WRITE(*,"(A,E14.4)")     "# HL      = ", HL
  WRITE(*,"(A,E14.4)")     "# HDL     = ", HDL
  WRITE(*,"(A,i14.4)")     "# NP     = ", NPR
  
!  call load_r(r_sponge)

!!$  OPEN(unit=91,file='energy.out') 
!!$  OPEN(unit=31, file='gulf.in',status='old') 
!!$  READ(31,*) 
!!$  READ(31,*) nsteps,dtf,ngraph,frame_int 
!!$  READ(31,*) 
!!$  READ(31,*) fnhhy 
!!$  READ(31,*) 
!!$  READ(31,*) dx,dy 
!!$  CLOSE(31) 
!!$  WRITE(6,*) 'nsteps= ',nsteps 
!!$
!!$  WRITE(6,*) 'parameters' 
!!$  WRITE(6,*) 'fnhhy=',fnhhy,'eps=',EPS,'delta=',delta,'qpr',qpr,    &
!!$       &     'lambda=',lambda,'beta=',beta,'dtime=',dtf,'WL',WL           
!!$  WRITE(6,*) 'R0, UL, P1, HL, HL/DL', R0,UL,P1,HL,HDL 

