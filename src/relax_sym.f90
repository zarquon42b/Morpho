
subroutine displace_mesh_gauss(iomat,nmat,Wvb,nw,Pvb,np,D1,D2,sigma,gamma,oway,clIW,ncl,clIP)
  IMPLICIT NONE
  integer :: nw, np,i,nmat,ncl,clIW(nmat,ncl),clIP(nmat,ncl),tmpW(ncl),tmpP(ncl)
  real*8 :: matr(nmat,3),Wvb(nw,3),Pvb(np,3),D1(nw,3),D2(np,3)
  real*8 :: sigma,gamma,k,iomat(nmat,3)
  real*8 :: point(3)
  logical :: oway
  !Pvb=Pvb+D2
  matr = iomat
! $OMP PARALLEL DO private(i, tmpW, tmpP,point) shared(Pvb, Wvb, ncl, iomat,gamma,sigma,oway)
  do i = 1,nmat
     tmpW(:) = clIW(i,:)
     tmpP(:) = clIP(i,:)
     point = matr(i,:)
     call relax_pt(point,Wvb(tmpW,:),ncl,Pvb(tmpP,:),ncl,D1(tmpW,:),D2(tmpP,:),sigma,gamma,iomat(i,:),oway)
     !outmat(i,1:3) = outpoint  
  end do
! $OMP END PARALLEL DO 

end subroutine displace_mesh_gauss

subroutine relax_pt(point,Wvb,nw,Pvb,np,D1,D2,sigma,gamma,outpoint,oway)
   IMPLICIT NONE
  integer :: nw, np,i
  real*8 :: point(3),Wvb(nw,3),Pvb(np,3),D1(nw,3),D2(np,3)
  real*8 :: sigma,gamma,k,g1tmp,g2tmp,g1sum,g2sum
  real*8 :: outpoint(3),outpoint2(3),diffvec(3)
  logical :: oway

  outpoint = outpoint*0
  outpoint2=outpoint
  g1sum = 0d0
  g2sum = 0d0
  
  do i = 1,nw     
     call gaussian_dsmooth(sigma,point,Wvb(i,:),g1tmp)
     outpoint(1:3) = outpoint(1:3) + g1tmp*D1(i,:)
     g1sum = g1sum + g1tmp
  end do
  outpoint = outpoint/g1sum
  
  if (oway .eqv. .false.) then
     do i = 1,np
        call gaussian_dsmooth(sigma,point,Pvb(i,:),g2tmp)
        outpoint2 = outpoint2 + g2tmp*D2(i,:)
        g2sum = g2sum + g2tmp
     end do
     outpoint2=outpoint2/g2sum
     outpoint = (outpoint-outpoint2)/gamma
  else
   outpoint = outpoint/gamma
  end if
 

end subroutine relax_pt

subroutine gaussian_dsmooth(sigma,point,tarpoint,out)

real*8 :: point(3), tarpoint(3),k,sigma,out,diffvec(3)

diffvec=point-tarpoint
k = dot_product(diffvec,diffvec)
out = exp(-k/sigma)

end subroutine gaussian_dsmooth

