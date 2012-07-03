subroutine barycenter(VB,nvb,IT,nit,bary)

  integer ::nit, IT(3,nit),i
  real*8 :: bary(nit,3),VB(3,nvb)

do i = 1,nit
   bary(i,1:3) = (VB(1:3,IT(1,i)) + VB(1:3,IT(3,i))+ VB(1:3,IT(2,i)))/3
end do 
end subroutine barycenter

subroutine face_zero(IT,nit,out)

  integer ::nit, IT(3,nit),i
  integer :: out(nit)

do i = 1,nit
   out(i) = product(IT(:,i))
end do 
end subroutine face_zero
