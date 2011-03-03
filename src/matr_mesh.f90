SUBROUTINE matr_mesh(matr,nmat,VB,nvb,IT,nit,dif,fptr,outmatr,regionv)

IMPLICIT NONE
integer ::nit, IT(3,nit),ittmp(3),i,nvb,nmat,fptr(nmat),ptrtmp,regionv(nmat)
integer:: region
real*8 :: point(3),clost(3),VB(3,nvb),normals(3,3),dif(nmat)
real*8 :: DAT(nit,12),matr(nmat,3),diff,outmatr(nmat,3)
 
 call updateSearch(VB,nvb,IT,nit,DAT)
 
 do i = 1,nmat

    point(:) = matr(i,:)
    
    call pt_upmesh(point,DAT,nit,clost,diff,ptrtmp,region)
    outmatr(i,:) = clost(:)
    dif(i) = diff
    fptr(i) = ptrtmp
    regionv(i) = region
    
 end do

END SUBROUTINE matr_mesh
