
subroutine displace_mesh_gauss(iomat,nm,Wvb,nw,Pvb,np,D1,D2,sigma,gam,ow,clIW,ncl,clIP,tol,rt0,rt1,rc)
!! displace a mesh according to a displacement vector field
  IMPLICIT NONE
!interface
!   subroutine relax_pt(pt,Wvb,nw,Pvb,np,D1,D2,sigma,gam,outpt,ow,dif1,dif2)
!     integer :: nw, np,i
!     real*8 :: pt(3),Wvb(nw,3),Pvb(np,3),D1(nw,3),D2(np,3)
!     real*8 :: sigma,gam,k,g1tmp,g2tmp,g1sum,g2sum
!     real*8 :: outpt(3),outpt2(3),diffvec(3),tmpd
!     logical :: ow, check
!     real*8, optional :: dif1(nw),dif2(np)
!   end subroutine relax_pt
!end interface

  integer :: nw, np,i,nm,ncl,clIW(nm,ncl),clIP(nm,ncl),tmpW(ncl),tmpP(ncl)
  real*8 :: matr(nm,3),Wvb(nw,3),Pvb(np,3),D1(nw,3),D2(np,3)
  real*8 :: sigma,gam,k,iomat(nm,3)
  real*8 :: pt(3),rt0(nw),rt1(np),rc
  real*8 :: dif1(nw),dif2(np)
  logical :: ow, check
  real*8 , intent(in) :: tol
  check = .false.
  matr = iomat
  where(dif1.ne.1)
     dif1 = 1 
  end where
  where(dif2.ne.1)
     dif2 = 1 
  end where
  
!! check if a threshold for displacement vector lengths is set
  if (tol > 0 .or. rc > 0) then
     
     check = .true.
     if (tol > 0) then !! calculate displacement lengths
        do i = 1,nw
           dif1(i) = sum(D1(i,:)**2)
        end do
        where (dif1 > tol) 
           dif1 = 0d0
        elsewhere
           dif1 = 1d0
        end where
        do i=1,np
              dif2(i) = sum(D2(i,:)**2)
           end do
           where (dif2 > tol) 
           dif2 = 0d0
        elsewhere
           dif2 = 1d0
        end where
     end if
     if (rc > 0) then !check for diverging normal angles
        where(rt0 > rc)
           dif1=0
        end where
        where(rt1 > rc)
           dif2=0
        end where
     end if
  end if

  do i = 1,nm
     tmpW(:) = clIW(i,:)
     tmpP(:) = clIP(i,:)
     pt = matr(i,:)
      
       
        call relax_pt(pt,Wvb(tmpW,:),ncl,Pvb(tmpP,:),ncl,D1(tmpW,:),D2(tmpP,:),sigma,gam,iomat(i,:),ow,dif1(tmpW),dif2(tmpP))
     
     !outmat(i,1:3) = outpt  
  end do

  
end subroutine displace_mesh_gauss
   !! calculate displacement weights according to moshfeghi, 1994
subroutine relax_pt(pt,Wvb,nw,Pvb,np,D1,D2,sigma,gam,outpt,ow,dif1,dif2)
  IMPLICIT NONE
  integer :: nw, np,i,tmpd1(nw),tmpd2(np)
  real*8 :: pt(3),Wvb(nw,3),Pvb(np,3),D1(nw,3),D2(np,3)
  real*8 :: sigma,gam,k,g1tmp,g2tmp,g1sum,g2sum
  real*8 :: outpt(3),outpt2(3),diffvec(3),tmpd
  logical :: ow, check
  real*8 :: dif1(nw),dif2(np)
  check = .false.
  
  !! initialize variables
  outpt = outpt*0
  outpt2=outpt
  g1sum = 0d0
  g2sum = 0d0
  !$OMP PARALLEL DO shared(sigma,Wvb,pt) private(i,g1tmp)
   !OMP& REDUCTION(+:g1sum,outpt) 
 
  do i = 1,nw  
     call gaussian_dsmooth(sigma,pt,Wvb(i,:),g1tmp)
     outpt(1:3) = outpt(1:3) + g1tmp*D1(i,:)*dif1(i)
     g1sum = g1sum + g1tmp*dif1(i)
  end do
  !$OMP END PARALLEL DO
  if (g1sum > 0) then
     outpt = outpt/g1sum
  end if
  
  if (ow .eqv. .false.) then
     do i = 1,np
        call gaussian_dsmooth(sigma,pt,Pvb(i,:),g2tmp)
        outpt2 = outpt2 + g2tmp*D2(i,:)*dif2(i)
        g2sum = g2sum + g2tmp*dif2(i)
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

