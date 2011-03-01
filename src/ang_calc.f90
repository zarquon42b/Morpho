SUBROUTINE angcal(xv,nx,yv,ny,rho)

 IMPLICIT NONE
	integer ::nx,ny
	real :: xv(nx),yv(ny),rho,lx,ly

lx = sqrt(dot_product(xv(:),xv(:)))
ly = sqrt(dot_product(yv(:),yv(:)))
 xv(:) = xv(:)/lx
 yv(:) = yv/ly
rho = (acos((sum((xv(:)-yv(:))**2)-2)))/(-2)


end subroutine angcal
