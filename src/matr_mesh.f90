SUBROUTINE matr_mesh(matr,nmat,VB,nvb,IT,nit,dif,fptr,outmatr,regionv,VBnormals,sign,outnorm)
  
  !! run brute force search of closest triangle on a mesh for each point given in 
  !! a matrix (matr)
  !! output: dif= distance; fptr=index of hit face, outmatr= matrix with closest points
  !! regionv= region of the triangle that was hit
  
  IMPLICIT NONE
  integer ::nit, IT(3,nit),ittmp(3),i,nvb,nmat,fptr(nmat),ptrtmp,regionv(nmat)
  integer :: region
  real*8 :: clost(3),VB(3,nvb),normals(3,3),dif(nmat),VBnormals(3,nvb),tmpnorm(3),tmpdiff(3)
  real*8 :: DAT(nit,13),diff,outmatr(nmat,3)
  real*8 ::  matr(nmat,3),outnorm(3,nmat)
  real*8 :: point(3),signo
  logical :: sign
  call updateSearch(VB,nvb,IT,nit,DAT)

! $OMP PARALLEL DO private(i, clost,tmpnorm,region,ptrtmp,diff) shared(outmatr,fptr,dif,regionv,outnorm,DAT)

  do i = 1,nmat
     
     point(:) = matr(i,1:3)
     
     call pt_upmesh(point,DAT,nit,clost,diff,ptrtmp,region)
     outmatr(i,:) = clost(:)
     dif(i) = diff
     fptr(i) = ptrtmp
     regionv(i) = region
     tmpnorm = VBnormals(:,IT(1,fptr(i)))+VBnormals(:,IT(2,fptr(i)))+VBnormals(:,IT(3,fptr(i)))
     tmpnorm = tmpnorm/sqrt(dot_product(tmpnorm,tmpnorm))
     outnorm(:,i) = tmpnorm 
     
     if (sign .eqv. .TRUE.) then
        
        tmpdiff = clost - point
        signo = dot_product(tmpdiff,tmpnorm)
        if ( signo < 0) then
           dif(i) = -dif(i)
        end if
        
        
     end if
  end do
! $OMP END PARALLEL DO 

  
END SUBROUTINE matr_mesh

SUBROUTINE matr_meshKD(matr,nmat,VB,nvb,IT,nit,clostInd,k,dif,fptr,outmatr,regionv,VBnormals,sign,outnorm)
  
  
  IMPLICIT NONE
  integer ::nit,k, IT(3,nit),ittmp(3),nvb, clostInd(nmat,k),nmat,fptr(nmat),ptrtmp,regionv(nmat),i
  integer :: region
  real*8 :: clost(3),VB(3,nvb),normals(3,3),dif(nmat)
  real*8 :: DAT(nit,13),diff,outmatr(nmat,3)
  real*8 :: matr(nmat,3),outnorm(3,nmat)
  real*8 :: point(3),signo,VBnormals(3,nvb),tmpnorm(3),tmpdiff(3)
  logical :: sign
  real*8 :: tmpdat(k,13)
  call updateSearch(VB,nvb,IT,nit,DAT)
! $OMP PARALLEL DO private(i, clost,tmpnorm,region,ptrtmp,diff) shared(outmatr,fptr,dif,regionv,outnorm,clostInd,DAT)
  
do i = 1,nmat
     
     point(:) = matr(i,1:3)
     tmpdat = DAT(clostInd(i,:),:)
     call pt_upmesh(point,tmpdat,k,clost,diff,ptrtmp,region)
     outmatr(i,:) = clost(:)
     dif(i) = diff
     fptr(i) = clostInd(i,ptrtmp)
     regionv(i) = region
     tmpnorm = VBnormals(:,IT(1,fptr(i)))+VBnormals(:,IT(2,fptr(i)))+VBnormals(:,IT(3,fptr(i)))
     tmpnorm = tmpnorm/sqrt(dot_product(tmpnorm,tmpnorm))
     outnorm(:,i) = tmpnorm
     
     if (sign .eqv. .TRUE.) then
        
        tmpdiff = clost - point
        signo = dot_product(tmpdiff,tmpnorm)
        if ( signo < 0) then
           dif(i) = -dif(i)
        end if
        
     end if
     
  end do
! $OMP END PARALLEL DO 

END SUBROUTINE matr_meshKD
