MODULE header

 integer, parameter :: rc_kind = selected_real_kind(12) 
 REAL(kind=rc_kind) :: toto3000
#include "cppdefs.h"

#include "size.h"

  integer :: NPR,ini_particle_time,parti_file_num

  LOGICAL, PARAMETER :: rect = .TRUE., periodicew = .TRUE.
!*** S0 and T0 are the mean salinity and temp. sreal= S0+s,Treal=T0+T  
  REAL(kind=rc_kind),  PARAMETER :: S0=35.7d0, T0=15.d0 ,R0=1027.d0
!*** EPS is the Rosby number
!*** OMEGA is the angular velocity of rotation of the earth
  REAL(kind=rc_kind),  PARAMETER :: EPS= 0.1d0, AL=1.d7, FPAR=1.d-4,                     &
                        PI=3.14159265358979323846, OMEGA=7.272d-5,           &
                        gpr= 0.981,apr=0.6371,dztop=1.d-3,z0= 0.2d0,zm= 0.1d0
!*** delta is the ratio of horizontal to vertical length scale = L/D   
!*** delinv= 1/delta
  REAL(kind=rc_kind) ::  delta,delinv
!*** L, DL are the characteristic length and depth scales.  
!*** lexp is a constant in the equations denoting the length scale     
!*** LEN= 10**(4+lexp),  0<=lexp<=1        
  REAL(kind=rc_kind), PARAMETER :: LEN= 1.d5, DL=1.d3
!*** dtime : time step
!*** dtf   : dimensionless time step
!*** fnhhy : 0 for hydrostatic, 1 for nonhydrostatic
  REAL(kind=rc_kind)  :: dtf,fnhhy,dx,dy,dtime
  integer :: pickup_step
!  REAL(kind=rc_kind)  :: dtf,fnhhy,dx,dy,dtime
!*** sigrelease is the isopycnal of tracer release
!*** dztop is the depth of the top layer non-dim by DL
!*** pfac is the grid stretching factor in z, used in findzall and sigma
  REAL(kind=rc_kind)  :: phi0deg
  REAL(kind=rc_kind)  :: yfront,dyfront,sclwidth,tightness,mldepth,total_depth
  REAL(kind=rc_kind)  :: stressmax,pfac,sigrelease(ntr)
  INTEGER :: conv(0:NI+1,0:NJ+1,0:NK+1)
  INTEGER :: con100(0:NI+1,0:NJ+1,0:NK+1)
  INTEGER :: ktest(3)
!*** qpr is the ratio of non-hyd to hyd pressure = Q/P
!*** qprinv = 1/qpr
!*** lambda is the ration of hor length scale to earth's rad = L/A 
!*** kappah is the implcitness parameter in the Crank-Nicolson scheme f
  REAL(kind=rc_kind)               :: sbackgrnd, beta,P1,qpr,lambda, &
        &                 UL,WL,TL,HL,HDL,DLinv,kappah,kaphinv

  REAL(kind=rc_kind) ::  Kx_TS, Ky_TS_bdy,Ky_TS_int, Kx_m, Ky_m
  REAL(kind=rc_kind) ::  r_sponge(0:NJ+1)

  REAL(kind=rc_kind), dimension(3) :: mldn2,mldn2init,zbtop,zbbot,frictop,diatop,     &
        &                 zbtopinit,zbbotinit,advecpv,friction,           &
        &                 diabatic,diabot,fricbot
 
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1,nconsume) :: consump 
  REAL(kind=rc_kind), target, dimension(ntr,0:NI+1,0:NJ+1, 0:NK+1,0:1)  :: Tr
  REAL(kind=rc_kind), pointer, dimension(:,:,:)  :: Tr_p
  REAL(kind=rc_kind), dimension(ntr,0:NI+1,0:NJ+1, 0:NK+1    )  :: wtr
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1,0:1)  :: u,v,w,s,T
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK,  3  )  :: gqk
  REAL(kind=rc_kind), dimension(    0:NI,    NJ,     NK,  2  )  :: gi,gqi
  REAL(kind=rc_kind), dimension(      NI,  0:NJ,     NK,  2  )  :: gj,gqj
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1,-1:NK+1)      :: zf
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1)      :: zc,wx,wy,wz,p,strain,shear,Jac,rho,&
        &                                              vor,pv,freqN2,si,sj,sk,cx,cy,cz,rp,T_ref
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK  )      :: wt,wzk,skfc
  REAL(kind=rc_kind), dimension(      NI,    NJ,   0:NK  )      :: czf,Kz,wf
  REAL(kind=rc_kind), dimension(           0:NJ+1,   NK  )      :: ueast,uwest
  REAL(kind=rc_kind), dimension(      NI,  0:NJ,     NK  )      :: hyn,sjfc,cyf,gj3,gqj3,Jjfc,vf
  REAL(kind=rc_kind), dimension(    0:NI,    NJ  ,   NK  )      :: hxn,sifc,cxf,gi3,gqi3,Jifc,uf
  REAL(kind=rc_kind), dimension(      NI,    NJ,     NK  )      :: uvis,vvis,wvis,fricu,fricv, &
        &                                              fricw,fricb,rhoadv,rhoprev
  REAL(kind=rc_kind), dimension(             NJ,     NK  )      :: ufbce,ufbcw,trinit,divreyn, &
        &                                              divmean,dcdt,prod
  REAL(kind=rc_kind), dimension(      NI,            NK  )      :: vfbcn,vfbcs,vnorth,vsouth,ssouth
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1,   2   )      :: gradhn
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1        )      :: ux,uy,vx,vy,ffc,bbc,oldh,h,hdt,D, &
        &                                              J2d,Ddx,Ddy,g11,g22,g12
  REAL(kind=rc_kind), dimension(      NI  ,0:NJ          )      :: bbj,ffj
  REAL(kind=rc_kind), dimension(    0:NI  ,  NJ          )      :: ffi,bbi
  REAL(kind=rc_kind), dimension(      NI,    NJ          )      :: wfbcb
  REAL(kind=rc_kind), dimension(           0:NJ+1        )      :: yc,latrad
  REAL(kind=rc_kind), dimension(    0:NI+1               )      :: xc
  REAL(kind=rc_kind), dimension(             NJ          )      :: stressx

!rpgrads
  REAL(kind=rc_kind) :: drpx(NI,NJ,NK),drpy(NI,NJ,NK),                                  &
            grpifc(0:NI,NJ,NK),grpjfc(NI,0:NJ,NK)

  CHARACTER(LEN=101) :: dirout

  ! added 20120302
 INTEGER :: ksurf, NP,it
 REAL(kind=rc_kind) :: hmean, drho, hlfac
 REAL(kind=rc_kind) ::  dum  
  

  ! added 20120302




 ! added 20120302
 ! netcdf output

  integer :: iddatfile, idigit, idudx, idudy, idudz, idvbysection, idvby, idvbz, idvcon, idvcy, idvbx
  integer :: idvcz, idvc, idvdivfeddy, idvdivfreyn, idvdx, idvdz, idvd, idvfb, idvh, idvn2bar
  integer :: idvn2, idvnsq100m, idvnsq30m, idvpe,idvpsiv,idvpsiw,idvpv, idvp,idvrhbar,idvrho,idvrnk
  integer :: idvstrain,idvstress,idvstr,idvs,idvtbar,idvtemp,idvtim,idvtr,idvt,idvu,idvvb,idvvc,idvvor,idvv
  integer :: idvwb,idvwc,idvwpv,idvw,idvy,idvzsave,idvz,idwdx,idwdy,idwdz,iimday,ipos





END MODULE header

