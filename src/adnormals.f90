SUBROUTINE adnormals(VB,M,IT,N,O,normals)

   IMPLICIT NONE
	integer :: M,i,j,t,N,IT(3,N),O
	real*8 :: VB(3,M),l,ntmp(1,1:3),tmp0(1,1:3),tmp1(1,1:3),normals(4,M),co
	data ntmp(1,1:3) / 0,0,0 /
	do i = 1,N
	tmp0(1,1:3) = VB(1:3,IT(1,i))-VB(1:3,IT(3,i))
	tmp1(1,1:3) = VB(1:3,IT(2,i))-VB(1:3,IT(1,i))
	!ntmp(1,1:3) = tmp0(1,1:3)
		call crossp(tmp0,tmp1,ntmp)
		call veclen(ntmp,l)
	if (l==0) then	
		ntmp(1,1:3) = (/0 , 0 , 0/)
		
		else
		ntmp = ntmp/l
	
		end if 	
		
		do j = 1,3
		co =sum(normals(1:3,IT(j,i))*ntmp(1,1:3))
		 if (co<0) then
		  normals(1:3,IT(j,i)) = normals(1:3,IT(j,i))-ntmp(1,1:3)
		  !call normalize(normals(1:3,IT(j,i)),l)
		  normals(4,IT(j,i)) = 1
		  
		 else 
		  normals(1:3,IT(j,i)) = normals(1:3,IT(j,i))+ntmp(1,1:3)
		   normals(1:3,IT(j,i)) = normals(1:3,IT(j,i))
		  !call normalize(normals(1:3,IT(j,i)),l)
		  normals(4,IT(j,i)) = 1
		 end if
		end do
 	end do
contains
 subroutine crossp (x,y,z)
	real*8 :: x(1,1:3),y(1,1:3),z(1,1:3)
	
	z(1,1) = x(1,2)*y(1,3)-x(1,3)*y(1,2)
	z(1,2) = x(1,3)*y(1,1)-x(1,1)*y(1,3)
	z(1,3) = x(1,1)*y(1,2)-x(1,2)*y(1,1)
 end subroutine crossp
 subroutine veclen (x,l)
	real*8 :: x(1,1:3),l
	l=sqrt(sum(x(1,1:3)**2))
	
 end subroutine veclen
subroutine normalize (x,l)
	real*8 :: x(1,1:3),l
	call veclen(x,l)	
	x(1,1:3) = x(1,1:3)/l
	
 end subroutine normalize
END SUBROUTINE adnormals


