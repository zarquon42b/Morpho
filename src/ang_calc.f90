SUBROUTINE angcal(xv,nx,yv,ny,rho,circle)

 IMPLICIT NONE
 integer :: nx,ny
integer:: circle
 real*8 :: xv(3),yv(3),rho,lx,ly, pi
! pi= ACOS(0.0)
!PI=4.D0*DATAN(1.D0)
 PI = 3.141592653589793239
 
 lx = sqrt(dot_product(xv(:),xv(:)))
 ly = sqrt(dot_product(yv(:),yv(:)))
!print *,lx
if (lx > 0) then 
 xv(:) = xv(:)/lx
end if
if (ly > 0) then
  yv(:) = yv(:)/ly
 end if
 rho = acos((sum((xv-yv)**2)-2)/(-2))

if (circle .eq. 0) then

   if (rho > pi/2) then
      rho = pi -  rho
   end if
end if
   
 !print *, xv

END SUBROUTINE angcal

subroutine angcheck(norm1,n,norm2,out,circle)
IMPLICIT NONE
integer :: n,i
integer,optional :: circle
real*8 :: norm1(3,n), norm2(3,n),out(n)
if (.not. present(circle)) then
   circle = 1
end if
do i = 1,n

 call angcal(norm1(:,i),3,norm2(:,i),3,out(i),circle)

end do

end subroutine angcheck

