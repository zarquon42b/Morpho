subroutine angcomp(angvec,ref,mat,i,j)

  integer :: i,j,k
  real*8 :: angvec(i),mat(i,j),ref(j),rho

  do k = 1,i
     call angcal(ref,j,mat(k,1:j),j,rho)
     angvec(k) = rho
  end do
  
  


end subroutine angcomp
