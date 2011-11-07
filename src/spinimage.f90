subroutine spinimage(Sp,ij,i,ab,imax,jmax,bs)
  implicit none
  integer :: ij(i,2),imax,jmax,k,i
  real*8 :: Sp(imax,jmax),ab(i,2),bs

do k = 1,i
   Sp(ij(k,1),ij(k,2)+1) = Sp(ij(k,1),ij(k,2)+1)+ab(k,1)*(bs-ab(k,2))
   Sp(ij(k,1),ij(k,2)) = Sp(ij(k,1),ij(k,2))+(bs-ab(k,1))*(bs-ab(k,2))        
   Sp(ij(k,1)+1,ij(k,2)+1) = Sp(ij(k,1)+1,ij(k,2)+1)+ab(k,1)*ab(k,2)
   Sp(ij(k,1)+1,ij(k,2)) = Sp(ij(k,1)+1,ij(k,2))+ab(k,2)*(bs-ab(k,1))
   
  

end do
end subroutine spinimage
