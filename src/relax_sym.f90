
subroutine displace_mesh_gauss(matr,nmat,Wvb,nw,Pvb,np,D1,D2,sigma,gamma,outmat)
   IMPLICIT NONE
  integer :: nw, np,i,nmat
  real*8 :: matr(nmat,3),Wvb(nw,3),Pvb(np,3),D1(nw,3),D2(np,3)
  real*8 :: sigma,gamma,k,outmat(nmat,3)
  real*8 :: outpoint(3),point(3),S0(nmat,3)
 
 
  do i = 1,nmat
     outpoint = outpoint*0
     point = matr(i,1:3)
     call relax_pt(point,Wvb,nw,Pvb,np,D1,D2,sigma,gamma,outpoint)
     outmat(i,1:3) = outpoint  
  end do


end subroutine displace_mesh_gauss

subroutine relax_pt(point,Wvb,nw,Pvb,np,D1,D2,sigma,gamma,outpoint)
   IMPLICIT NONE
  integer :: nw, np,i
  real*8 :: point(3),Wvb(nw,3),Pvb(np,3),D1(nw,3),D2(np,3)
  real*8 :: sigma,gamma,k,g1tmp,g2tmp,g1sum,g2sum
  real*8 :: outpoint(3),outpoint2(3)
  !outpoint = outpoint*0
  g1sum = 0
  g2sum = 0

  do i = 1,nw
     
     call gaussian_one(sigma,point,Wvb(i,:),g1tmp)
     outpoint = outpoint + g1tmp*D1(i,:)
     g1sum = g1sum + g1tmp
  end do
  outpoint = outpoint/g1sum
  
  do i = 1,np
     call gaussian_two(sigma,point,Pvb(i,:),D2(i,:),g2tmp)
     !call gaussian_one(sigma,point,Pvb(i,:),g2tmp)

     outpoint2 = outpoint2 + g2tmp*D2(i,:)
     g2sum = g2sum + g2tmp
  end do
  outpoint2=outpoint2/g2sum
  outpoint = (outpoint-outpoint2)/gamma
 

end subroutine relax_pt

subroutine gaussian_one(sigma,point,tarpoint,out)


real*8 :: point(3), tarpoint(3),k,sigma,out,diffvec(3)

diffvec=point-tarpoint
k = dot_product(diffvec,diffvec)
out = exp(-k/sigma)

end subroutine gaussian_one

subroutine gaussian_two(sigma,point,tarpoint,displace,out)
  IMPLICIT NONE
real*8 :: point(3), tarpoint(3),displace(3),k,sigma,out,diffvec(3)
diffvec=point-tarpoint-displace
k = dot_product(diffvec,diffvec)
out = exp(-k/sigma)

end subroutine gaussian_two
