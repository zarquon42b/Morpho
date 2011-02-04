SUBROUTINE tpsfx(A,M,B,N,Bh,x,coefs,tot)

   IMPLICIT NONE
	integer :: M,i,t,N
	real*8 :: A(M,3),tmp(1,1:3),B(N,3),coefs(3,M+4),x(M),naf,tot(N,3),coefna(3,M),Bh(N,4)
	coefna = coefs(1:3,1:M)
	!Bh(N,1) = 1
	!Bh(N,2:4) = B(N,1:3)
	!data x /5*0/
	do t= 1,N   	!calculate nonaffine part
	
		do i = 1,M
  		tmp(1,1:3) = A(i,1:3) - B(t,1:3)
		x(i) = -sqrt(sum(tmp(1,1:3) **2))
		!x(i) = matmul(coefs,)!x(i)*coefs(i)
		!naf = sum(x)
		end do
	tmp(1,1:3) = matmul(coefna,x)
	
	tot(t,1:3) = matmul(coefs(1:3,(M+1):(M+4)),Bh(t,1:4))+tmp(1,1:3)
	
	end do
	
END SUBROUTINE tpsfx


