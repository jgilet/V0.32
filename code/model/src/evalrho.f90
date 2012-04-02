subroutine evalrho(rhonew,n) 
!     ---------------------------------------------                     
  USE header,only:s,T,NI,NJ,NK,rc_kind
!  s is salinity, T is temp, rhonew is pot density
  implicit none
  REAL(kind=rc_kind) :: rhonew(0:NI+1,0:NJ+1,0:NK+1),potdens
  integer i,j,k,n 
  do k=0,NK+1 
     do j=0,NJ+1 
        do i=0,NI+1 
           rhonew(i,j,k)= potdens(s(i,j,k,n),T(i,j,k,n))
           !rhonew(i,j,k)= rho(i,j,k)
        end do
     end do
  end do
  return 
END subroutine evalrho
