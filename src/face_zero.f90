subroutine face_zero(IT,nit,out)

  integer :: nit, IT(3,nit),i
  integer :: out(nit)

do i = 1,nit
   out(i) = product(IT(:,i))
end do 
end subroutine face_zero
