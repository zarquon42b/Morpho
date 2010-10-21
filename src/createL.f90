SUBROUTINE createL(K,k1,M,n)

   	IMPLICIT NONE
	integer :: 	k1,i,j,n
	real*8 :: 	K(k1,k1), M(k1,n)
		do i = 1,k1
			do j = 1,k1
			K(i,j) = sqrt(sum((M(i,1:n)-M(j,1:n))**2))
			end do
		end do
	
	
END SUBROUTINE createL
