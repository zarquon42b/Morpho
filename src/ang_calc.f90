SUBROUTINE angcal(xv,nx,yv,ny,rho)

 IMPLICIT NONE
 integer :: nx,ny
 real*8 :: xv(3),yv(3),rho,lx,ly

 lx = sqrt(dot_product(xv(:),xv(:)))
 ly = sqrt(dot_product(yv(:),yv(:)))
!print *,lx
 xv(:) = xv(:)/lx
 yv(:) = yv(:)/ly
 
 rho = acos((sum((xv-yv)**2)-2)/(-2))
 !print *, xv

END SUBROUTINE angcal
