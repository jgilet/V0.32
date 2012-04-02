subroutine energy(step) 
!---------------------------------------------------                    
  USE header
                                                                       
      integer :: step 
      REAL(kind=rc_kind) ::  const,enkin,enkin2
                                                                       
      const= EPS*EPS*delta*delta 

      enkin = SUM( Jac(1:NI,1:NJ,1:NK) *           &
                  ( u(1:NI,1:NJ,1:NK,0)**2 +         &
                    v(1:NI,1:NJ,1:NK,0)**2 +         &
                    const* w(1:NI,1:NJ,1:NK,0)**2  ))
              
!      enkin2= 0.d0 
!      do k=1,NK 
!         do j=1,NJ 
!            do i=1,NI 
!               enkin2= Jac(i,j,k)*(u(i,j,k,0)*u(i,j,k,0) +v(i,j,k,0)*    &
!     &              v(i,j,k,0) + const*w(i,j,k,0)*w(i,j,k,0))           &
!     &              +enkin                                              
!            end do 
!         end do 
!      end do 

      print*, "#total kinetic energy = ", step, enkin
!      print*, "#total kinetic energy = ", step, enkin2
      return 
END subroutine energy                                       
