subroutine meanMat(A,n,m,out)
integer :: n,m,i
real*8::A(n,m), out(m)

do i=1,m
   out(i) = sum(A(:,i))/n
end do


end subroutine meanMat
