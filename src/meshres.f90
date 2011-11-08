subroutine meshres(VB,nvb,IT,nit,res)

  integer :: nvb,nit,i,j,IT(3,nit)
  real*8 :: VB(3,nvb),res,tmp(3)

  do i = 1,nit
     tmp = VB(1:3,IT(1,i))-VB(1:3,IT(2,i))
     res= res+sqrt(dot_product(tmp,tmp))
     tmp = VB(1:3,IT(1,i))-VB(1:3,IT(3,i))
     res= res+sqrt(dot_product(tmp,tmp))
     tmp = VB(1:3,IT(2,i))-VB(1:3,IT(3,i))
     res= res+sqrt(dot_product(tmp,tmp))
  end do
  res = res/(nit*3)

end subroutine meshres
