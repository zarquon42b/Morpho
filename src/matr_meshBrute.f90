SUBROUTINE matr_meshBrute(matr,nmat,VB,nvb,IT,nit,dif,fptr,outmatr,regionv,VBnormals,sign,outnorm,mm)
  
  !! run brute force search of closest triangle on a mesh for each point given in 
  !! a matrix (matr)
  !! output: dif= distance; fptr=index of hit face, outmatr= matrix with closest points
  !! regionv= region of the triangle that was hit
  
  IMPLICIT NONE
  integer ::nit, IT(3,nit),ittmp(3),i,nvb,nmat,fptr(nmat),ptrtmp,regionv(nmat)
  integer :: region,mm, j
  real*8 :: clost(3),VB(3,nvb),normals(3,3),dif(nmat),VBnormals(3,nvb),tmpnorm(3),tmpdiff(3)
  real*8 :: DAT(nit,13),diff,outmatr(nmat,3)
  real*8 ::  matr(nmat,3),outnorm(3,nmat)
  real*8 :: point(3),signo,nlen, weight(3)
  logical :: sign
  call updateSearch(VB,nvb,IT,nit,DAT)

! $OMP PARALLEL DO private(i, clost,tmpnorm,region,ptrtmp,diff) shared(outmatr,fptr,dif,regionv,outnorm,DAT)

  do i = 1,nmat
     
     point(:) = matr(i,1:3)
     
     call pt_upmesh(point,DAT,nit,clost,diff,ptrtmp,region,mm)
     outmatr(i,:) = clost(:)
     dif(i) = diff
     fptr(i) = ptrtmp
     regionv(i) = region
     do j = 1, 3
        tmpdiff = clost-VB(:,IT(j,fptr(i)))
        weight(j) = sqrt(dot_product(tmpdiff,tmpdiff))
        if (weight(j) .eq. 0) then
           weight(j) = 1e12
        else
           weight(j) = 1/weight(j)
        end if
     end do
     weight = weight/sum(weight)
     tmpnorm = weight(1)*VBnormals(:,IT(1,fptr(i)))+ weight(2)*VBnormals(:,IT(2,fptr(i)))+ weight(3)*VBnormals(:,IT(3,fptr(i)))
     nlen=sqrt(dot_product(tmpnorm,tmpnorm))
     if (nlen > 0) then
        tmpnorm = tmpnorm/nlen
     end if
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

  
END SUBROUTINE matr_meshBrute
