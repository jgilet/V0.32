subroutine stprofile 
  !     --------------------                                              
!  USE header
  ! Initializes s,T from Andrey's section gridded for model grid
  ! s_init.dat and T_init.dat  written in ascii - NK+2 x NJx2
!  implicit none
!  integer n,i,j,k

  !     --------------------                                              
  USE header, ONLY: NI,NJ,NK,rho,s,T,pi,zc,DL,LEN,yc,dy,dx,zf
  implicit none
  INTEGER n,i,j,k,Nu,k0
  REAL*8 :: A,B,G,C,potdens,x,y,z(0:NK),tmp,yfront,mldepth,width,slfac,anom
  INTEGER,PARAMETER :: seed = 86456  
  CALL RANDOM_SEED()

  ! assign the density or temperature and salinity by either analytic functions or 
  ! any particular hydrographic section.

  ! t=A/B exp(A*z)erf(B*y) + G*(1+z)
  ! z=-1:0, y = -0.5,0.5, 

  ! ===========================
  ! cases 
  ! 0: jet
  ! 1: idealized surface jet
  ! 2: krushio
  SELECT CASE (1)
  CASE(1)
     A = 4.0d0
     B = 40d3 !jet has 2B width
     C = 300d0
     G = 19d0
     !Nu =NK-9
     mldepth = 50d0
     s(:,:,:,0) = 34d0

     DO k = 0,NK
        z(k) = (zc(1,1,k) )*DL
        anom = tanh((z(k))/C)*(1-tanh((z(k))/C)**2)
        print*, anom
        T(:,:,k,0) = 30d0 - 10d0 * tanh(z(k)/C)**2
        DO j = 0, NJ+1
              y = -(dble(j)-dble(NJ+1)/2d0)*dy
           DO i = 0, NI+1
              x = -(dble(i)-dble(NI+1)/2d0)*dx
        T(i,j,k,0) = T(i,j,k,0) - anom* exp(-(y**2+(x-z(k)/DL*100e3)**2)/B**2)
              rho(i,j,k) = potdens(s(i,j,k,0),T(i,j,k,0))
     ENDDO
     ENDDO
     ENDDO

     T(:,0,:,0)= T(:,1,:,0)
     T(:,NJ+1,:,0)= T(:,NJ,:,0)
     T(:,:,NK+1,:)=T(:,:,NK,:)
     s(:,0,:,0)= s(:,1,:,0)
     s(:,NJ+1,:,0)= s(:,NJ,:,0)
     s(:,:,NK+1,:)=s(:,:,NK,:)
     rho(:,0,:)= rho(:,1,:)
     rho(:,NJ+1,:)= rho(:,NJ,:)
     rho(:,:,NK+1)=rho(:,:,NK)
!     s(:,:,:,0) = rho 
!     DO k = 1, NK
!        rho(:,:,k) = 0.25d0*(rho(:,:,k-1)+2d0*rho(:,:,k)+rho(:,:,k+1))
!        s(:,:,k,0) = 0.25d0*(s(:,:,k-1,0)+2d0*s(:,:,k,0)+s(:,:,k+1,0))
!        T(:,:,k,0) = 0.25d0*(T(:,:,k-1,0)+2d0*T(:,:,k,0)+T(:,:,k+1,0))
!     ENDDO
     call save3d(NI+2,NJ+2,NK+2,rho(:,:,:),'initial-rho.bin')
  END SELECT

END SUBROUTINE stprofile
