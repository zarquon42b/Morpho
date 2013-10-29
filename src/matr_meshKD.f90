SUBROUTINE matr_meshKD(matr,nmat,VB,nvb,IT,nit,clostInd,k,dist,fptr,regionv,VBnormals,sign,outnorm,mm,baryc,baryn)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  closest point search on a triangle surface mesh                !!!
!!!  based on a search of the k-closted triangles determined by an  !!!
!!!  intitial kd-tree search across faces' barycenters              !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IMPLICIT NONE
  integer ::nit,k, IT(3,nit),ittmp(3),nvb, clostInd(nmat,k),nmat,fptr(nmat),ptrtmp,regionv(nmat),i
  integer :: region,mm,j, baryn
  real*8 :: clost(3),VB(3,nvb),normals(3,3),dist(nmat)
  real*8 :: DAT(nit,13),diff, weight(3)
  real*8 :: matr(nmat,3),outnorm(3,nmat)
  real*8 :: point(3),signo,VBnormals(3,nvb),tmpnorm(3),tmpdiff(3)
  logical :: sign
  real*8 :: tmpdat(k,13),nlen, baryc(3, baryn)
  ! prepare search structures
  call updateSearch(VB,nvb,IT,nit,DAT)
! $OMP PARALLEL DO private(i, clost,tmpnorm,region,ptrtmp,diff) shared(outmatr,fptr,dist,regionv,outnorm,clostInd,DAT)
  
do i = 1,nmat
     
     point(:) = matr(i,1:3)
     tmpdat = DAT(clostInd(i,:),:)
     ! run search on k-subset determined by KD-tree search
     call pt_upmesh(point,tmpdat,k,clost,diff,ptrtmp,region,mm)
     matr(i,:) = clost(:)
     dist(i) = diff
     fptr(i) = clostInd(i,ptrtmp)
     regionv(i) = region
     ! get normal weights from surrounding vertices - inverse distance
     do j = 1, 3
        tmpdiff = clost-VB(:,IT(j,fptr(i)))
        weight(j) = sqrt(dot_product(tmpdiff,tmpdiff))
        if (weight(j) .eq. 0) then
           weight(j) = 1e12
        else
           weight(j) = 1/weight(j)
        end if
     end do
     ! apply weighting to normals
     tmpnorm = weight(1)*VBnormals(:,IT(1,fptr(i)))+ weight(2)*VBnormals(:,IT(2,fptr(i)))+ weight(3)*VBnormals(:,IT(3,fptr(i)))
     ! normalize normals
     nlen=sqrt(dot_product(tmpnorm,tmpnorm))
     if (nlen > 0) then
        tmpnorm = tmpnorm/nlen
     end if
     outnorm(:,i) = tmpnorm
     ! get signed distance
     if (sign .eqv. .TRUE.) then
        tmpdiff = clost - point
        signo = dot_product(tmpdiff,tmpnorm)
        if ( signo < 0) then
           dist(i) = -dist(i)
        end if
     end if
     ! get barycentric coordinates
     if (baryn .eq. nmat) then 
        call barycoo(baryc(:,i),clost,fptr(i) , VB, nvb, IT,nit)
     end if
     
     
  end do
  ! $OMP END PARALLEL DO 
 contains 
   SUBROUTINE barycoo(barytmp,point,fptrtmp,VB, nvb, IT,nit)
     integer :: nit, nvb
     integer :: IT(3,nit),fptrtmp
     real*8 :: barytmp(3), point(3), d00, d01, d11, d20, d21
     real*8 :: VB(3,nvb), v0(3), v1(3), v2(3), denom
     v0=VB(:,IT(2,fptrtmp))-VB(:,IT(1,fptrtmp))
     v1=VB(:,IT(3,fptrtmp))-VB(:,IT(1,fptrtmp))
     v2=point - VB(:,IT(1,fptrtmp))
     d00 = dot_product(v0,v0)
     d01 = dot_product(v0,v1)
     d11 = dot_product(v1,v1)
     d20 = dot_product(v2,v0)
     d21 = dot_product(v2,v1)
     denom = d00*d11 -d01**2
     barytmp(2) =  (d11 * d20 - d01 * d21) / denom
     barytmp(3) = (d00 * d21 - d01 * d20) / denom
     barytmp(1) = 1 - barytmp(2) - barytmp(3)
   end SUBROUTINE barycoo
END SUBROUTINE matr_meshKD
