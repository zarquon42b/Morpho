SUBROUTINE angcal(xv,nx,yv,ny,rho)

 IMPLICIT NONE
 integer :: nx,ny
 real*8 :: xv(3),yv(3),rho,lx,ly

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
 !print *, xv

END SUBROUTINE angcal

subroutine angcheck(norm1,n,norm2,out)
IMPLICIT NONE
integer :: n,i 
real*8 :: norm1(3,n), norm2(3,n),out(n)

do i = 1,n

 call angcal(norm1(:,i),3,norm2(:,i),3,out(i))

end do

end subroutine angcheck

