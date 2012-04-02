subroutine advection(m,n,dtime,step) 
#include "cppdefs.h"
  USE header, only : NI,NJ,NK
  !     ---------------------------------------------                     
  !     convecn(m,n,dtime) means use level m and write s,T to level n     
  !     computes Cx,Cy,Cz at the cell centers.                            
  !     uflx,vflx,wflx,sflx,Tflx are the fluxes defined at                
  !     cell faces (using QUICK).                                         
  !                                                                       
  integer i,j,k,n,m,step,it 
  !     ntr is the number of tracers in variable T (Defined in header.f)  
  real*8 ::absvel,left,right,ctr,upos,uneg 
  real*8 ::dtime,uflx(0:NI,0:NJ,0:NK),                      &
       vflx(0:NI,0:NJ,0:NK),wflx(0:NI,0:NJ,0:NK),                   &
       sflx(0:NI,0:NJ,0:NK),Trflx(ntr,0:NI,0:NJ,0:NK),              &
       Tflx(0:NI,0:NJ,0:NK),uTx(0:NI+1,0:NJ+1,0:NK+1),              &
       usx(0:NI+1,0:NJ+1,0:NK+1),uTrx(ntr,0:NI+1,0:NJ+1,0:NK+1),    &
       uux(0:NI+1,0:NJ+1,0:NK+1),uvx(0:NI+1,0:NJ+1,0:NK+1),         &
       uwx(0:NI+1,0:NJ+1,0:NK+1),Tdif(NI,NJ,NK),                    &
       sdif(NI,NJ,NK),Trdif(ntr,NI,NJ,NK),udif(NI,NJ,NK),           &
       vdif(NI,NJ,NK),wdif(NI,NJ,NK),dtJ                            
  real*8 ::dTrx(ntr,0:NI+1),dsx(0:NI+1),dux(0:NI+1),        &
       dTx(0:NI+1),dvx(0:NI+1),dwx(0:NI+1),dTry(ntr,0:NJ),          &
       dTy(0:NJ),dsy(0:NJ),duy(0:NJ),dvy(0:NJ),dwy(0:NJ),           &
       dTrz(ntr,0:NK+1),dsz(0:NK+1),dTz(0:NK+1),                    &
       duz(0:NK+1),dvz(0:NK+1),dwz(0:NK+1)                          
  real*8 ::s_restore(0:NJ+1),alph_restore(0:NJ+1) 
  real*8 ::2d0,fc 

  real*8, dimension(0:NI+1,0:NJ+1,0:NK+1) :: dTx
  real*8, dimension(0:NI,0:NJ,0:NK) :: dTx
  real*8, dimension(1:NI,1:NJ,1:NK) :: uf,vf,wf,upos,uneg,left,right,ctr
                                                                        
      m=0
      dTx(0:NI,1:NJ,1:NK) = T(1:NI+1,1:NJ,1:NK,m) - & 
                            T(0:NI  ,1:NJ,1:NK,m)
      dTx(NI+1,1:NJ,1:NK) = dTx(1,1:NJ,1:NK)

      absvel= dabs(uf) 
      upos= 0.5d0*(uf +absvel) 
      uneg= 0.5d0*(uf -absvel) 
                                                                        
      left  = 0.125d0* (dTx(1:NI,1:NJ,1:NK)   - dTx(0:NI-1,1:NJ,1:NK)) 
      right = 0.125d0* (dTx(2:NI+1,1:NJ,1:NK) - dTx(1:NI,1:NJ,1:NK)) 
      ctr   = 0.5d0  * (T(2:NI+1,1:NJ,1:NK,m) + T(1:NI,1:NJ,1:NK,m)) 
      Tflx(1:NI,1:NJ,1:NK) = uf*ctr - upos*left -uneg*right 
      Tflx(0,1:NJ,1:NK)  = Tflx(NI,1:NJ,1:NK)
      uTx(1:NI,1:NJ,1:NK) = Tflx(1:NI,1:NJ,1:NK) - Tflx(0:NI-1,1:NJ,1:NK)  
                                                                        
!!      do k=1,NK 
!!         do j=1,NJ 
!     for periodic-ew boundaries                                        
!     dTx,dsx,.... at faces...                                          
!          i=0 
!           dTx(i)= T(1,j,k,m) - T(NI,j,k,m) 
!           do i=1,NI-1 
!              dTx(i)= T(i+1,j,k,m) - T(i,j,k,m) 
!           end do 
!           i= NI 
!           dTx(i)= T(1,j,k,m) - T(i,j,k,m) 
!     specially for periodic-ew boundaries dsx,dTx... from 0:NI+1       
!           i= NI+1 
!           dTx(i)= T(2,j,k,m) - T(1,j,k,m) 
!     periodic-ew boundaries  : do i=1,NI  (instead of i=1,NI-1)        
!           T(NI+1,j,k,m)= T(1,j,k,m) 
                                                                        
!!            do i=1,NI 
!!               goto 108 
!!               if (uf(i,j,k).gt.0.d0) then 
!!                                                                        
!!                  left= 0.125d0*(dTx(i) -dTx(i-1)) 
!!                  ctr= 0.5d0* (T(i+1,j,k,m) +T(i,j,k,m)) 
!!                  fc= ctr -left 
!!                  call ultim(T(i-1,j,k,m),T(i,j,k,m),                   &
!!     &                 T(i+1,j,k,m),dTx(i),dTx(i-1),2d0,fc)         
!!                  Tflx(i,j,k)= uf(i,j,k)*fc 
!!               else 
!!
!!                  right= 0.125d0*(dTx(i+1) -dTx(i)) 
!!                  ctr= 0.5d0* (T(i+1,j,k,m) +T(i,j,k,m)) 
!!                  fc= ctr - right 
!!                  call ultim(T(i+2,j,k,m),T(i+1,j,k,m),                 &
!!     &                 T(i,j,k,m),dTx(i),dTx(i+1),2d0,fc)           
!!                  Tflx(i,j,k)= uf(i,j,k)*fc 
!!               end if 
                                                                        
!     THIS PART was without the limiter ultim in x-dir                  
!==               goto 109                                              
!!  108          continue 
!!               absvel= dabs(uf(i,j,k)) 
!!               upos= 0.5d0*(uf(i,j,k) +absvel) 
!!               uneg= 0.5d0*(uf(i,j,k) -absvel) 
!!                                                                        
!!               left= 0.125d0*(dTx(i) -dTx(i-1)) 
!!               right= 0.125d0*(dTx(i+1) -dTx(i)) 
!!               ctr= 0.5d0* (T(i+1,j,k,m) +T(i,j,k,m)) 
!!               Tflx(i,j,k) = uf(i,j,k)*ctr - upos*left -uneg*right 
!!                                                                        
!!  109       end do 
!     periodic-ew boundaries                                            
 !!           Tflx(0,j,k)= Tflx(NI,j,k) 
 !!        end do 
 !!     end do 


!     do k=1,NK 
!        do j=1,NJ 
!           do i=1,NI 
!              uTx(i,j,k)=  Tflx(i,j,k) -Tflx(i-1,j,k) 
!           end do 
!        end do 
!     end do 
                                                                        
!     y-direction                                                       
!     -----------                                                       
      do k=1,NK 
         do i=1,NI 
!     dTy,dsy,.... at faces...                                          
            do j=0,NJ,NJ 
               dTy(j)= 0.d0 
            end do 
            do j=1,NJ-1 
               dTy(j)= T(i,j+1,k,m) - T(i,j,k,m) 
            end do 
            T(i,0,k,m)= T(i,1,k,m) 
                                                                        
            T(i,NJ+1,k,m)= T(i,NJ,k,m) 
                                                                        
            do j=1,NJ-1 
                                                                        
               if (vf(i,j,k).ge.0.d0) then 
                                                                        

                  left= 0.125d0*(dTy(j) -dTy(j-1)) 
                  ctr= 0.5d0* (T(i,j+1,k,m) +T(i,j,k,m)) 
                  fc= ctr -left 
                  call ultim(T(i,j-1,k,m),T(i,j,k,m),                   &
     &                 T(i,j+1,k,m),dTy(j),dTy(j-1),2d0,fc)         
                  Tflx(i,j,k)= vf(i,j,k)*fc 
                                                                        
               else 
                                                                        
                  right= 0.125d0*(dTy(j+1) -dTy(j)) 
                  ctr= 0.5d0* (T(i,j+1,k,m) +T(i,j,k,m)) 
                  fc= ctr - right 
                  call ultim(T(i,j+2,k,m),T(i,j+1,k,m),                 &
     &                 T(i,j,k,m),dTy(j),dTy(j+1),2d0,fc)           
                  Tflx(i,j,k)= vf(i,j,k)*fc 
                                                                        
               end if 
                                                                        
            end do 
!     Solid Boundaries                                                  
            do j=0,NJ,NJ 
               Tflx(i,j,k)= 0.d0 
            end do 
         end do 
      end do 
                                                                        
      do k=1,NK 
         do i=1,NI 
            do j=1,NJ 
               uTx(i,j,k)= (Tflx(i,j,k) -Tflx(i,j-1,k) ) + uTx(i,j,k) 
            end do 
         end do 
      end do 
                                                                        
!     z-direction                                                       
!     -----------                                                       
      do j=1,NJ 
         do i=1,NI 
!     dTz,dsz,.... at faces...                                          
!            do k=0,NK,NK                                               
            k=0 
            dTz(k)= 0.d0 
!                                                                       
            T(i,j,k,m)= T(i,j,k+1,m) 
!                                                                       
            do k=1,NK-1 
               dTz(k)= T(i,j,k+1,m) - T(i,j,k,m) 
            end do 
!     Top boundary - Zero gradient                                      
            do k=NK,NK+1 
               dTz(k)= 0.d0 
            end do 
!     Linear extrapolation                                              
            T(i,j,NK+1,m)= T(i,j,NK,m) 
                                                                        
            do k=1,NK 
                                                                        
               if (wf(i,j,k).ge.0.d0) then 
                                                                        
                  left= 0.125d0*(dTz(k) -dTz(k-1)) 
                  ctr= 0.5d0* (T(i,j,k+1,m) +T(i,j,k,m)) 
                  fc= ctr -left 
                  call ultim(T(i,j,k-1,m),T(i,j,k,m),                   &
     &                 T(i,j,k+1,m),dTz(k),dTz(k-1),2d0,fc)         
                  Tflx(i,j,k)= wf(i,j,k)*fc 
               else 
                                                                        
                  if (k.eq.NK) then 
                     Tflx(i,j,k)= wf(i,j,k)*0.5d0*                       &
     &                    (T(i,j,k+1,m) +T(i,j,k,m))                    
                     goto 209 
                  endif 
                                                                        

                  right= 0.125d0*(dTz(k+1) -dTz(k)) 
                  ctr= 0.5d0* (T(i,j,k+1,m) +T(i,j,k,m)) 
                  fc= ctr - right 
                  call ultim(T(i,j,k+2,m),T(i,j,k+1,m),                 &
     &                 T(i,j,k,m),dTz(k),dTz(k+1),2d0,fc)           
                  Tflx(i,j,k)= wf(i,j,k)*fc 
                                                                        
               end if 
  209       end do 
!     Solid Boundary                                                    
            k=0 
            Tflx(i,j,k)= 0.d0 
         end do 
      end do 
                                                                        
      call diffusion(sdif,Trdif,udif,vdif,m,step) 
                                                                        
!=      call biharmonic(udif,vdif,m,step)                               
!1/25/05 moved     call viscous(udif,vdif,wdif,m)                       
!                                                                       
      do 250 k=1,NK 
         do 251 j=1,NJ 
            do 252 i=1,NI 
               uTx(i,j,k)= uTx(i,j,k) + (Tflx(i,j,k) -Tflx(i,j,k-1) )   &
     &              - Tdif(i,j,k)                                       
!     put back sdif to calc. advec of rho, non-dim by Tl                
               rhoadv(i,j,k)= (usx(i,j,k) +sdif(i,j,k))/Jac(i,j,k) 
  252       continue 
  251    continue 
  250 continue 
      call viscous(udif,vdif,wdif,sdif,Tdif,Trdif,m) 
      do k=1,NK 
         do j=1,NJ 
            do i=1,NI 
               uTx(i,j,k)= uTx(i,j,k)-Tdif(i,j,k) 
            end do 
         end do 
      end do 
!                                                                       
!                                                                       
!     *$*  ASSERT CONCURRENT CALL                                       
      do 300 j=1,NJ 
         do 301 i=1,NI 
            do 302 k=1,NK 
               dtJ= dtime/Jac(i,j,k) 
               T(i,j,k,n)= T(i,j,k,0) -dtJ*uTx(i,j,k) 
!==     &              - alph_restore(j)*(s(i,j,k,m)-s_restore(j))      
!     T_restore= 0.d0                                                   
  302       continue 
  301    continue 
  300 continue 
!                                                                       
!+      call remineralize(n) ( called from momentum)                    
!     Surface flux                                                      
!=      call surfaceflux(dtime,n)                                       
!                                                                       
      return 
      END                                           
                                                                        
                                                                        
      subroutine ultim(up,ct,dn,gdn,gup,2d0,fc) 
!     ----------------------------------------------                    
!     Implementation of the ULTIMATE limiter based on upstrm,dnstrm,ctr 
!     points and d(phi)/dx at the up- and dnstrm locations              
!     Returns the "ulitimate limited" face value of the variable        
!     that should be multiplied by uf,vf or wf to give the flux.        
                                                                        
      real*8 ::up,ct,dn,gup,gdn,2d0,fc,ref,del,adel,acurv 
                                                                        
      del= dn -up 
      adel= dabs(del) 
      acurv= dabs(gdn -gup) 
      if (acurv.ge.adel)  then 
         fc= ct 
         return 
      else 
         ref= up + (ct -up)*2d0 
         if (del.gt.0) then 
            fc= dmax1(ct,fc) 
            ref= dmin1(ref,dn) 
            fc= dmin1(fc,ref) 
         else 
            fc= dmin1(ct,fc) 
            ref= dmax1(ref,dn) 
            fc= dmax1(ref,fc) 
         endif 
      end if 
      return 
      END                                           
