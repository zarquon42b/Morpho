 subroutine crossp (x,y,z)
	real*8 :: x(1,1:3),y(1,1:3),z(1,1:3)
	
	z(1,1) = x(1,2)*y(1,3)-x(1,3)*y(1,2)
	z(1,2) = x(1,3)*y(1,1)-x(1,1)*y(1,3)
	z(1,3) = x(1,1)*y(1,2)-x(1,2)*y(1,1)
 end subroutine crossp
