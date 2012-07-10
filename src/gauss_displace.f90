
subroutine displace_mesh_gauss(iomat,nm,Wvb,nw,Pvb,np,D1,D2,sigma,gam,ow,clIW,ncl,clIP,tol)
!! displace a mesh according to a displacement vector field
  IMPLICIT NONE
interface
   subroutine relax_pt(pt,Wvb,nw,Pvb,np,D1,D2,sigma,gam,outpt,ow,dif1,dif2)
     integer :: nw, np,i
     real*8 :: pt(3),Wvb(nw,3),Pvb(np,3),D1(nw,3),D2(np,3)
     real*8 :: sigma,gam,k,g1tmp,g2tmp,g1sum,g2sum
     real*8 :: outpt(3),outpt2(3),diffvec(3),tmpd
     logical :: ow, check
     real*8, optional :: dif1(nw),dif2(np)
   end subroutine relax_pt
end interface

  integer :: nw, np,i,nm,ncl,clIW(nm,ncl),clIP(nm,ncl),tmpW(ncl),tmpP(ncl)
  real*8 :: matr(nm,3),Wvb(nw,3),Pvb(np,3),D1(nw,3),D2(np,3)
  real*8 :: sigma,gam,k,iomat(nm,3)
  real*8 :: pt(3)
  real*8, allocatable:: dif1(:),dif2(:)
  logical :: ow, check
  real*8 , intent(in) :: tol
  check = .false.
  matr = iomat

  !! check if a threshold for displacement vector lengths is set
  if (tol > 0 ) then
     !! calculate displacement lengths
     allocate(dif1(nw))
     allocate(dif2(np))
     check = .true.
     do i = 1,nw
        dif1(i) = sum(D1(i,:)**2)
        if (dif1(i) > tol) then
           dif1(i) = 0d0
        end if
     end do
     do i=1,np
        dif2(i) = sum(D2(i,:)**2)
        if (dif2(i) > tol) then
           dif2(i) = 0d0
        end if
     end do
  end if

  ! $OMP PARALLEL private(i, tmpW, tmpP,pt) shared(Pvb, Wvb, ncl, iomat,gam,sigma,ow)
  !OMP DO
  do i = 1,nm
     tmpW(:) = clIW(i,:)
     tmpP(:) = clIP(i,:)
     pt = matr(i,:)
     if (check .eqv. .true.) then   
       
        call relax_pt(pt,Wvb(tmpW,:),ncl,Pvb(tmpP,:),ncl,D1(tmpW,:),D2(tmpP,:),sigma,gam,iomat(i,:),ow,dif1(tmpW),dif2(tmpP))
     else
        call relax_pt(pt,Wvb(tmpW,:),ncl,Pvb(tmpP,:),ncl,D1(tmpW,:),D2(tmpP,:),sigma,gam,iomat(i,:),ow)
     end if
     !outmat(i,1:3) = outpt  
  end do
! OMP END DO
  ! $OMP END PARALLEL DO 
  
end subroutine displace_mesh_gauss
   !! calculate displacement weights according to moshfeghi, 1994
subroutine relax_pt(pt,Wvb,nw,Pvb,np,D1,D2,sigma,gam,outpt,ow,dif1,dif2)
  IMPLICIT NONE
  integer :: nw, np,i
  real*8 :: pt(3),Wvb(nw,3),Pvb(np,3),D1(nw,3),D2(np,3)
  real*8 :: sigma,gam,k,g1tmp,g2tmp,g1sum,g2sum
  real*8 :: outpt(3),outpt2(3),diffvec(3),tmpd
  logical :: ow, check
  real*8, optional :: dif1(nw),dif2(np)
  check = .false.
  
  if (present(dif1)) then
     check= .true.
  end if
  !! initialize variables
  outpt = outpt*0
  outpt2=outpt
  g1sum = 0d0
  g2sum = 0d0
  
  do i = 1,nw  
     tmpd = 1d0
     if (check .eqv. .true.) then
        tmpd = dif1(i)
     end if
     if (tmpd > 0) then
        call gaussian_dsmooth(sigma,pt,Wvb(i,:),g1tmp)
        outpt(1:3) = outpt(1:3) + g1tmp*D1(i,:)
        g1sum = g1sum + g1tmp
     end if
  end do
  if (g1sum > 0) then
     outpt = outpt/g1sum
  end if
  
  if (ow .eqv. .false.) then
     do i = 1,np
        tmpd = 1d0
        if (check .eqv. .true.) then
           tmpd = dif2(i)
        end if
        if (tmpd > 0) then
           call gaussian_dsmooth(sigma,pt,Pvb(i,:),g2tmp)
           outpt2 = outpt2 + g2tmp*D2(i,:)
           g2sum = g2sum + g2tmp
        end if
     end do
     if (g2sum > 0) then
        outpt2=outpt2/g2sum
     end if
     outpt = (outpt-outpt2)/gam
  else
     outpt = outpt/gam
  end if
end subroutine relax_pt

subroutine gaussian_dsmooth(sigma,pt,tarpt,out)
!! calculate specific weight according to radial basis function  
  real*8 :: pt(3), tarpt(3),k,sigma,out,diffvec(3)
  
  diffvec=pt-tarpt
  k = dot_product(diffvec,diffvec)
  out = exp(-k/sigma)
  
end subroutine gaussian_dsmooth

