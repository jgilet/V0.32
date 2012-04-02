PROGRAM main 
  !---------------------------------------------------                    
  !     test model
#include "cppdefs.h"

  USE header
  use relaxation
  use grids
#ifdef particle
  USE particles
#endif

  !======================================================================
  ! initionlization part, to be moved to a seperate file
#include "ini_param.f90"
! open(unit=91,file=TRIM(dirout)//'energy.out')
 call set_coef()
 call ini_grids()

#ifdef particle
  !======================================================================
  ! allocate partical array after NP is signed
  ALLOCATE(parti(NPR))
#endif
  !======================================================================
  ! define pointer, map the 5-dimensional Tr 
  ! to the 3-dimensional Tr_p for saving
  Tr_p => Tr(1,:,:,:,0)
!  allocate(pcorr(maxout))
  !======================================================================
  !      call potfunc                                                     
  !      call init(p,vfent,ufex)                                          

   PRINT*,'============='
   PRINT*,rc_kind
   PRINT*,'============='


  CALL init(pcorr) 
  !==      call stprofile (now called from init_tr.f)                     
  CALL meanh(NI,NJ,h,hmean) 
  WRITE(*,*) 'initial  hmean',hmean
  !-      call stinit    -- for variation with latitude                   
  CALL sigma           ! calculates metric terms in vertical for the moving part of the grid
  CALL staticsigma     ! calculates metric terms in vertical for the fixed part of the grid
  !-      call staticsigma (move the whole grid)                          
  !     important to call sigma before staticsigma and stprofile          
  !=      call stprofile                                                  
  !=      call bioinit                                                    

  CALL tracerinit()    !initializes tracer
  !      call setbcuf                                                     
  CALL hsave            ! saves initial free surface h 
  !     hsave must be called before geostroph                             
  CALL geostroph(0)     ! find geostrophically balanced velocities from T,s,h
  !     sigma needs to be called in momentum after each time advancement  
  !     check the divergence of the initial vel field                     

  CALL facediv(EPS,fdiv)          ! checks the divergence using face fluxes
  CALL cdiv(EPS,ctrdiv,0)         ! checks the divergence using face fluxes
  !  WRITE(*,*) ' fdiv  ',fdiv,' cdiv  ',ctrdiv
  WRITE(*,*) '# init called'
  !     write out initial values                                          
  CALL vort(0)                 ! calculates vorticity  -diagnostic
  PRINT*,"W ",MAXVAL(w)

  step= 0 
  CALL calcn2 
  CALL n2budget(step) 
  WRITE(*,*) "# to start momentum"
  PRINT*, "# dtime =",dtime
  tim= dtime
  !     BCs set just once in the beginning                                
  CALL setbc(step) 
  CALL correctbc 
  PRINT*,"W ",MAXVAL(w)

  CALL checks 
  PRINT*,"W ",MAXVAL(w),MINVAL(w)

!call w_pickup('op.pickup.'//stepchar(step)//'.bin')
!pickup files
!  if (pickup_step > 0) then
!call r_pickup(pickup_step)
!  end if
ksurf=NK;call writeksurf(frame_int,step,ksurf,h,consump,Tr,s,T,rho,u,v,w,p,vor,strain,freqN2,xc,yc,zc,DL,LEN,Jac,dtf*TL)
ksurf= INT(NK/3); call writeksurf(frame_int,step,ksurf,h,consump,Tr,s,T,rho,u,v,w,p,vor,strain,freqN2,xc,yc,zc,DL,LEN,Jac,dtf*TL)
PRINT*," step",step
  PRINT*,"W ",MAXVAL(w),MINVAL(w)
    CALL energy(step) 

  DO step = pickup_step+1,pickup_step+nsteps
#ifdef particle
     ! initialize particles
     if (step == ini_particle_time) then
        CALL ini_particles(step)
       !get the particle velocity for the first step
       CALL get_parti_vel(step)
     DO i = 1, NP
      parti(i)%i = parti(i)%i + dtf * parti(i)%u
      parti(i)%j = parti(i)%j + dtf * parti(i)%v
      parti(i)%k = parti(i)%k + dtf * parti(i)%w
      parti(i)%u0 = parti(i)%u
      parti(i)%v0 = parti(i)%v
      parti(i)%w0 = parti(i)%w
  ENDDO
   endif
#endif

   
     !convert integer 'step' into character 'stepchar' with zero padding on the left.
     !stepchar has the length of 10, which is supposely large enough for reasonable runs.
!     write(stepchar,'(I10.10)') step
!     call stepchar(step)
     PRINT*, 'steps='//stepchar(step)

    CALL momentum(pcorr,step) 
     !     calcn2 must be called before n2budget. Called at end of momentum  
    CALL n2budget(step) 
    !     compute the KE                                                    
    CALL energy(step) 
    CALL meanh(NI,NJ,h,hmean) 

#ifdef particle
     if (step>ini_particle_time) then
     !particles
     CALL get_parti_vel(step)
     CALL parti_forward
     call save_parti
     endif
#endif

#include "write_op.f90"
  ENDDO
! close(91)
END PROGRAM main
