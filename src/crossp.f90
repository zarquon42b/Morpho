 subroutine crossp (x,y,z)
	real*8 :: x(3),y(3),z(3)
	
	z(1) = x(2)*y(3)-x(3)*y(2)
	z(2) = x(3)*y(1)-x(1)*y(3)
	z(3) = x(1)*y(2)-x(2)*y(1)
 z = z/sqrt(dot_product(z,z))
 end subroutine crossp
