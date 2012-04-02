subroutine conadjust(n) 
  !----------------------------------------------------                   
  USE header
  !     T(i,j,k) is set to 0, so Treal is T0                              
  !     Performs convective adjustment                                    
  integer i,j,k,n,it,step 
  real*8 ::rh,rhtop,zt,zb,dz,zinv 
  real*8 ::dzupper 
  real*8 ::sreal,Treal,rhreal,pbar 
!      real*8 ::rho(0:NI+1,0:NJ+1,0:NK+1)                       
!                                                                       
!-      call evalrho(rho,n)                                             
                                                                        
      do j=1,NJ 
         do i=1,NI 
!                                                                       
            do k=NK-1,1,-1 
               conv(i,j,k)= 0 
                                                                        
               rhtop= s(i,j,k+1,n) 
!     top                                                               
               dzupper = zf(i,j,k+1) -zf(i,j,k) 
                                                                        
!     this level                                                        
               dz = zf(i,j,k) -zf(i,j,k-1) 
               rh =  s(i,j,k,n) 
                                                                        
               if (rh.lt.rhtop) then 
!     mix                                                               
                  conv(i,j,k)= 1 
                  zinv= 1.d0/(dzupper +dz) 
                  s(i,j,k,n)=(dzupper*s(i,j,k+1,n) +dz*s(i,j,k,n))*zinv 
                  s(i,j,k+1,n)= s(i,j,k,n) 
                  do it=1,ntr 
                     T(it,i,j,k,n)=(dzupper*T(it,i,j,k+1,n)             &
     &                    +dz*T(it,i,j,k,n))*zinv                       
                     T(it,i,j,k+1,n)= T(it,i,j,k,n) 
                  end do 
               end if 
                                                                        
            end do 
         end do 
      end do 
                                                                        
      if (mod((step-1),100).eq. 0) then 
         do k=0,NK+1 
            do j=0,NJ+1 
               do i=0,NI+1 
                  con100(i,j,k) = 0 
               end do 
            end do 
         end do 
      else 
         do k=0,NK+1 
            do j=0,NJ+1 
               do i=0,NI+1 
                  con100(i,j,k) = con100(i,j,k) + conv(i,j,k) 
               end do 
            end do 
         end do 
      end if 
!     re-evaluate density (if needed)                                   
!      call evalrho(rho,n)                                              
                                                                        
!-      call sTbc_periodicew(n)                                         
      return 
      END                                           
