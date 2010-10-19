SUBROUTINE tpsfx(A,M,B,N,x,coefs,tot)

   IMPLICIT NONE
	integer :: M,i,t,N
	real*8 :: A(M,3),tmp(1,1:3),B(N,3),coefs(M+4),x(M),naf,tot(N)
	
	!data x /5*0/
	do t= 1,N   	!calculate nonaffine part
	
		do i = 1,M
  		tmp(1,1:3) = A(i,1:3) - B(t,1:3)
		x(i) = sqrt(sum(tmp(1,1:3) **2))
		x(i) = x(i)*coefs(i)
		naf = sum(x)
		end do
	naf = sum(x)
	tot(t) = coefs(M+1)+coefs(M+2)*B(t,1)+coefs(M+3)*B(t,2)+coefs(M+4)*B(t,3)+naf	
	end do
END SUBROUTINE tpsfx


