subroutine spinmap(SO,VB,p,n,dimvb)
  implicit none
  integer :: i,dimvb
  real*8 :: SO(dimvb,2),VB(3,dimvb),p(3),x(3),n(3)

do i = 1,dimvb
   x = VB(1:3,i)
   x = x - p
   SO(i,1) = sqrt(sum(x**2)-dot_product(n,x)**2)
   SO(i,2) = dot_product(n,x)

end do
end subroutine spinmap
