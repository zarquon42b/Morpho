SUBROUTINE mesh_mesh(matr,nmat,VB,nvb,IT,nit,dif,fptr,outmatr,regionv,matnorm,tarnorm,rhotol,toldist)!(matr,matnorm,nmat,VB,tarnorm,nvb,IT,nit,dif,fptr,outmatr)

  IMPLICIT NONE
  integer ::nit, IT(3,nit),ittmp(3),i,nvb,nmat,fptr(nmat),ptrtmp,regionv(nmat)
  integer :: region,j
  real*8 :: clost(3),normal(3),dif(nmat),lnorm,rhotol,toldist,normtmp(3)
  real*8 :: VB(3,nvb),DAT(nit,12),diff,tarnorm(3,nvb),itnorm(nit,3)!target mesh
  real*8 :: matr(nmat,3),matnorm(nmat,3)!input matrix + normals
  real*8 :: point(3),outmatr(nmat,3)!output
 ! toldist = 1
 ! rhotol = 0.5
  !call updateSearch(VB,nvb,IT,nit,DAT)

  do i = 1,nit
     ittmp(:) = IT(:,i)
     normtmp(:) = tarnorm(:,ittmp(1)) + tarnorm(:,ittmp(2)) + tarnorm(:,ittmp(3))

     lnorm = sqrt(dot_product(normtmp,normtmp))
     normtmp = normtmp/lnorm
     itnorm(i,:) = normtmp(:)
     
    
  end do
! print *,itnorm(1,:)
!  do j = 1,nmat

 !    point(:) = matr(j,1:3)
 !    normal(:) = matnorm(j,1:3)
 !    call pt_upmesh(point,DAT,nit,clost,diff,ptrtmp,region)
     !call pt_angmesh(point,normal,DAT,nit,itnorm,clost,diff,ptrtmp,rhotol)
  !   outmatr(j,:) = clost(:)
  !   dif(j) = diff
  !   fptr(j) = ptrtmp
     !regionv(i) = region

  !end do
  call updateSearch(VB,nvb,IT,nit,DAT)
!   print *,nit
!   print *,nmat


   !iterate over matrix rows

   do i = 1,nmat

      point(:) = matr(i,1:3)
      normal(:) = matnorm(i,:)
      !call pt_upmesh(point,DAT,nit,clost,diff,ptrtmp,region)
      call pt_angmesh(point,normal,DAT,nit,itnorm,clost,diff,ptrtmp,rhotol,toldist)
      outmatr(i,:) = clost(:)
      dif(i) = diff
      fptr(i) = ptrtmp
      regionv(i) = region

   end do
 END SUBROUTINE mesh_mesh
