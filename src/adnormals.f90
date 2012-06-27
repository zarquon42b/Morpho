SUBROUTINE adnormals(VB,nvb,IT,nit,O,normals)

   IMPLICIT NONE
	integer :: nvb,i,j,t,nit,IT(3,nit),O
	real*8 :: VB(3,nvb),l,ntmp(1,1:3),tmp0(1,1:3),tmp1(1,1:3),normals(4,nvb),co
	data ntmp(1,1:3) / 0,0,0 /
	do i = 1,nit
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
		  normals(4,IT(j,i)) = normals(4,IT(j,i)) + 1
		  
		 else 
		  normals(1:3,IT(j,i)) = normals(1:3,IT(j,i))+ntmp(1,1:3)
		   normals(1:3,IT(j,i)) = normals(1:3,IT(j,i))
		  !call normalize(normals(1:3,IT(j,i)),l)
		  normals(4,IT(j,i)) = normals(4,IT(j,i)) + 1
		 end if
		end do
 	end do
	do i = 1,nvb
	normals(1:4,i) = normals(1:4,i)/normals(4,i)
	call veclen(normals,l)
	call normalize(normals(1:3,i),l)
	end do
contains
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


