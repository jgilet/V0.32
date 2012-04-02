SUBROUTINE init(pcorr)  
  !     subroutine init(p,vfent,ufex)                                    
  !     ------------------------------------------------                  
  !     For the curvilinear grid.                                         
  !     sigma levels are evenly spaced.  Hence wz is a function of        
  !     time but not of z.   Here u,v,w refer to the xi,eta,sigma direcito
  !     Those metric quantities that do not change with time are evaluated
  !     The rest are evaluated in "sigma.f" which will be called at every 
  !     step.                                                             
  !     Also it is absolutely essential that the integral of the flux     
  !     over the entrance section is equal to that over the exit section. 
  !     -No longer necessary with the free surface.                       
  !     This may require some special attention when the free surface is  
  !     moving at the exit. At the entrance vf is held fixed.             
  !                                                                       
  USE header 
  use grids

  !      REAL(kind=rc_kind) :: lat(0:NI+1,0:NJ+1),lon(0:NI+1,0:NJ+1)           
  !      REAL(kind=rc_kind) :: lat(0:NJ+1)                                     
  REAL(kind=rc_kind) ::   xdu(0:NI+1,0:NJ+1),ydu(0:NI+1,0:NJ+1),                       &
       &     xdv(0:NI+1,0:NJ+1),ydv(0:NI+1,0:NJ+1),C,fconst,              &
       &     phi0,cosphi0,sinphi0,sumuf,sumvf,dthet,dtheta,dphi,          &
       &     gin,c1,c2,temp,temp2,dnkinv,tmplon,tmplat,cjent,             &
       &     pcorr(maxout),dep,y 

  ! added 20120302
  REAL(kind=rc_kind) ::  dscl, rdscl, yarg,ex2y,thy,sechy2,ran3

  !                                                                       
  !     FPAR is the value used for scaling Coriolis parameters f and b    
  INTEGER i,j,k,jmid,iseed 

  !      open(unit=40, file='xyD.dat')                                    
  !      open(unit=45, file='du.metrics.dat')                             
  !      open(unit=50, file='dv.metrics.dat')                             
  !      open(unit=55, file='partialsD.dat')                              
  !      open(unit=110,file='ufopen.in')                                  
  WRITE(*,*) 'in init' 
  DLinv= 1.d0/DL 
  C= PI/180.d0 
  gin= 1.d0/gpr 
  dnkinv= 1.d0/DBLE(NK) 
  !                                                                       
  !     dx,dy are the grid dims (in m) read in in main.f                  
  !     phi0 is the central lat, dphi,dtheta grid spacing(angle)          
  !=      phi0deg= 25.d0                                                  
!  phi0deg=42.d0 
  !      phi0deg= 3.d0                                                    
  phi0= phi0deg*C 
  print*, "# phi0deg = ",phi0deg
  dphi = dy/(apr*AL) 
  print*, '# dphi = ',dphi, 'in deg ',dphi/C 
  !     The domain should be within the latitude limits  in stinit.f.     
  !     lat0 is the center latitude about which we linearize              
  jmid= NJ/2 
  print*,'# jmid=',jmid 
  DO  j=0,NJ+1 
     latrad(j)= phi0 + DBLE(j-jmid)*dphi 
     !latrad(j)= phi0                                               
  enddo
     !      do 10 j=0,NJ+1                                                   
     !         do 20 i=0,NI+1                                                
     !            read(40,*) lon(i,j),lat(i,j),D(i,j)                        
     !            lat(i,j)= lat(i,j)*C                                       
     !     rescale lat and lon so that the grid spacing is about 1/8 degree  
     !            lon(i,j)= lon(i,j)*C                                       
     ! 20      continue                                                      
     ! 10   continue                                                         
     !      close(40)                                                        
     !c      cosphi0= dCos(phi0)                                             
     !c      sinphi0= dSin(phi0)                                             
     !c      write(6,*) 'cosphi0=',cosphi0                                   
     !      lon0= 260.d0*C                                                   
     !      lat0= 20.d0*C                                                    
     !      dtheta= C/8.d0                                                   
     !      dphi= C/8.d0                                                     
     !      phi0= lat0 + 0.5d0*(dphi*dfloat(NJ))                             
     !c     for the rectilinear case, dx,dy,dz are constant                  
     !      dx= AL*apr*cosphi0*dtheta /LEN                                   
     !      dy= AL*apr*dphi/LEN                                              
     !      dz= 1000.d0/(dfloat(NK)*DL)                                      
     !      write(6,*) 'dx,dy,dz',dx,dy,dz                                   
     !      do 10 i=-1,NI+1                                                  
     !         do 10 j=-1,NJ+1                                               
     !            lat(i,j)= lat0 + dphi*dfloat(j)                            
     !            lon(i,j)= lon0 + dtheta*dfloat(i)                          
     ! 10   continue                                                         
     !     These values must be read in for the curvilinear grid.            
     !c      c1= apr*AL*cosphi0/LEN                                          
     !c      c2= apr*AL/LEN                                                  

!read dy information from dyM.data   
!     call load_dy(dyM)
! print*, dyM

     DO  j=0,NJ+1 
        DO  i=0,NI+1 
           xdu(i,j)= dx/LEN 
           ydu(i,j)= 0.d0 
           xdv(i,j)= 0.d0 
            ydv(i,j)= dy/LEN !- used for constant dy
!           ydv(i,j)= dyM(j)/LEN 
        ENDDO
     ENDDO
     !                                                                       
!---- the following block is used for varying dy -moved to mod_grids.f90 -J.Wang 21-Dec-2011
!     yc(0) = -0.5*dyM(0)*1.d-3 
!     DO j=1,NJ+1 
!        yc(j)= yc(j-1) + dyM(j)*1.d-3 
!     END DO
!---- the following block is used for constant dy
!     yc(0) = -0.5*dy*1.d-3 
!     DO j=1,NJ+1 
!        yc(j)= yc(j-1) + dy*1.d-3 
!     END DO
     xc(0) = -0.5*dx*1.d-3 
     DO i=1,NI+1 
        xc(i)= xc(i-1) + dx*1.d-3 
     END DO
     DO  j=0,NJ+1 
        DO  i=0,NI+1 
           J2d(i,j)= xdu(i,j)*ydv(i,j) -xdv(i,j)*ydu(i,j) 
           ux(i,j)= ydv(i,j)/J2d(i,j) 
           vx(i,j)= -ydu(i,j)/J2d(i,j) 
           uy(i,j)= -xdv(i,j)/J2d(i,j) 
           vy(i,j)= xdu(i,j)/J2d(i,j) 
           g11(i,j)= ux(i,j)*ux(i,j) +uy(i,j)*uy(i,j) 
           g12(i,j)= ux(i,j)*vx(i,j) +uy(i,j)*vy(i,j) 
           g22(i,j)= vx(i,j)*vx(i,j) +vy(i,j)*vy(i,j) 
        ENDDO
     ENDDO
     !                                                                       
     !     D(i,j) is -ve and  non-dim by DL                                  
     dep = total_depth
     print*, "total depth",dep
     D(:,:)= -dep*DLinv
     !     Compute partial  dD/dx,dD/dy                                      
     CALL smooth 
     WRITE(6,*) 'smooth called' 

     ! 101  call findzall                                                    
     !101  CONTINUE
     !      close(55)                                                        
     !c      open (unit= 31, file='xyD.dat')                                 
     !c      do j=0,NJ+1                                                     
     !c         tmplat= phi0deg + dble(j-jmid)*dthet                         
     !c         yc(j)= tmplat                                                
     !c         do i=0,NI+1                                                  
     !c            tmplon= 300.d0 +dthet*dble(i)                             
     !c            write(31,*) tmplon,tmplat,D(i,j)                          
     !c            if (j.eq.0) xc(i)= tmplon                                 
     !c         end do                                                       
     !c      end do                                                          
     !c      close(31)                                                       
     !     In smooth we smooth D and then use central differencing to evaluat
     !     Ddx and Ddy. The values obtained are much lower than those obtaine
     !     by smoothing Ddx, Ddy evaluated from the bspline - probably bec   
     !     we had a bug with a factor of 3 which is now corrected. (values   
     !     differ even after the correction and we use the numerial values.) 
     !c      call smooth                                                     
     !*    flat bottom                                                       
     ! 999  do j=0,NJ+1                                                      
     !         do i=0,NI+1                                                   
     !            Ddx(i,j)= 0.d0                                             
     !            Ddy(i,j)= 0.d0                                             
     !            D(i,j)= 1.d0                                               
     !            write(150,*) Ddx(i,j)                                      
     !            write(250,*) Ddy(i,j)                                      
     !            write(350,*) D(i,j)                                        
     !         end do                                                        
     !      end do                                                           
     !      stop                                                             
     !                                                                       
     u(:,:,:,0) = 0d0
     v(:,:,:,0) = 0d0
     w(:,:,:,0) = 0d0
     s(:,:,:,0) = 0d0
     T(:,:,:,0) = 0d0
     Tr(:,:,:,:,0) = 0d0
     conv(:,:,:)    = 0
     con100(:,:,:)  = 0

     pcorr(:)= 0.d0 
     print*, "###"

     uvis(:,:,:) = 0.d0
     vvis(:,:,:) = 0.d0
     wvis(:,:,:) = 0.d0
     !                                                                       
     !     specify the initial fields                                        
     !     read in the cell-centered pressure and u,v velocities             
     h = 0d0
     uf = 0d0
     vf = 0d0
     wf = 0d0
     !     fill in the vfbc arrays                                           

     ufbce(:,:)= uf(NI,:,:) 
     ufbcw(:,:)= uf(0,:,:) 


     vfbcn(:,:)= vf(:,NJ,:) 
     vfbcs(:,:)= vf(:,0,:) 

     !     at the sea bed, wf=0                                              

     wfbcb(:,:) = 0.d0 
     !     Evaluate the Coriolis parameter f' and fill the arrays ffi,ffj,ffc
     !     f= f0+By= 2Omega(Sin(phi0) +(phi-phi0)Cos(phi0))                  
     !     b= b0+By= 2Omega(Cos(phi0) -(phi-phi0)Sin(phi0))                  
     fconst= 2.d0*OMEGA/FPAR 

    
     DO j =1, NJ
        ffi(:,j)= fconst*sin(latrad(j)) 
        bbi(:,j)= fnhhy*fconst*cos(latrad(j)) 
     ENDDO
     DO j =0, NJ
        ffj(:,j)= fconst*sin(0.5d0*(latrad(j+1)+latrad(j))) 
        bbj(:,j)= fnhhy*fconst*cos(0.5d0*(latrad(j+1)+latrad(j))) 
     ENDDO
     DO j=0, NJ+1
        ffc(:,j)= fconst*sin(latrad(j)) 
        bbc(:,j)= fnhhy*fconst*cos(latrad(j)) 
     ENDDO
     !                                                                       
     !     Initialize s,T                                                    
     CALL findzall   ! finds the vertical grid
     print*,zc(10,10,:)
     print*, 'findzall'
!     CALL ini_rho()
     call stprofile
     call evalrho(rho,0) 
    print*, 'after ini_rho()'
     CALL inith           ! initialize free surface h
     print*, 'hmean',hmean
     CALL findzall        ! find z again
     !     write out z-grid                                                  
     !     -----------------                                                 
     OPEN (unit=60,file=TRIM(dirout)//'zgrid.out')
     WRITE(60,*) '# vertical grid' 
     DO k=0,NK+1 
        WRITE(60,*) k,zc(10,10,k)*1000. 
     END DO
     WRITE(60,*) '# face values' 
     DO k=-1,NK+1 
        WRITE(60,*) k,zf(10,10,k)*1000. 
     END DO
     CLOSE(60) 
     !     -----------                                                       


     advecpv(:)= 0.d0 
     friction(:)= 0.d0 
     diabatic(:)= 0.d0 
!if (pickup_step == 0) then
     !     Perturb the front boundary if the model does not start from pickup                                     
     !     for a quicker onset of instability                                
     !     ----------------------------------------------------              
!     iseed= 44294 
!     dum = ran3(iseed) 
!     !                                                                       
!     DO i=0,NI+1 
!        perturb =  1.d-3*ran3(iseed)
!        DO j=1,NJ 
!           h(i,j)= h(i,j) +perturb 
!        END DO
!     END DO
!endif

     RETURN 

     drho= rho(NI/2,5*NJ/6,NK)-rho(NI/2,NJ/6,NK) 
     print*, '# drho = ',drho 
     hlfac = 0.5*drho*100./DL/HL 
     !     100 is the depth of the front : upper ocean  depth                

     !     ADD a PERTURBATION to the width of the front                      
     iseed= 44294 
     dum = ran3(iseed) 

     !cc      hlfac=.1d0/HL  (set above)                                     
     DO i=0,NI+1 
        ! perturb                                                               
        dscl=0.1d0*(0.5*(yc(NJ)+yc(NJ+1))-0.5*(yc(0)+yc(1))) 
        dscl=(1.d0 + 1.d-4*ran3(iseed))*dscl 
        rdscl=1.d0/dscl 
        dscl=dscl*1.d-2 
        IF (i.EQ.10) OPEN(65) 
        !         write(6,*) 'dscl',dscl                                        
        DO j=0,NJ+1 
           !     ulfac=gpr*hlfac/(dscl*ffc(0,j))                                   
           yarg=yc(j)-  0.5*(yc(nj/2)+yc(nj/2+1)) 
           ex2y=EXP(-2.d0*rdscl*ABS(yarg)) 
           thy=(1.d0-ex2y)/(1.d0+ex2y) 
           IF(yarg.LT.0.d0) thy=-thy 
           sechy2=1.d0-thy*thy 

           hdt(i,j)= 0.d0 
           h(i,j)= -hlfac*thy 
        END DO
        IF (i.EQ.10) WRITE(65,*) j,h(10,j),u(10,j,15,0) 
     END DO
     CLOSE(65) 

     RETURN 

     hlfac=.1d0/HL 
     dscl=0.1d0*(0.5*(yc(NJ)+yc(NJ+1))-0.5*(yc(0)+yc(1))) 
     !cmy      dscl=0.25d0*(0.5*(yc(NJ)+yc(NJ+1))-0.5*(yc(0)+yc(1)))         
     rdscl=1.d0/dscl 
     dscl=dscl*1.d-2 
     OPEN(65) 
     DO j=0,NJ+1 
        !         ulfac=gpr*hlfac/(dscl*ffc(0,j))                               
        yarg=yc(j)-  0.5*(yc(nj/2)+yc(nj/2+1)) 
        ex2y=EXP(-2.d0*rdscl*ABS(yarg)) 
        thy=(1.d0-ex2y)/(1.d0+ex2y) 
        IF(yarg.LT.0.d0) thy=-thy 
        sechy2=1.d0-thy*thy 
        DO i=0,NI+1 
           hdt(i,j)= 0.d0 
           h(i,j)= -hlfac*thy 
           !            do k=0,nk+1                                                
           !               u(i,j,k,0)=ulfac*sechy2                                 
           !            enddo                                                      
        END DO
        WRITE(65,*) j,h(10,j),u(10,j,15,0) 
     END DO
     CLOSE(65) 

     RETURN 

     !                                                                       
     !     call evalrho(0)                                                   
     !      write(6,*) 'in initsq, evalrho called'                           
     !      write(6,*) 'writing out layer depths'                            
     !      do k=0,NK                                                        
     !         temp= findz(dble(k),dztop,D(10,10),0.d0,NK)                   
     !         write(6,*) 'k,z',k,temp                                       
     !      end do                                                           
     WRITE(6,*) 'to leave init' 
     RETURN 
     ! 
   END SUBROUTINE init
