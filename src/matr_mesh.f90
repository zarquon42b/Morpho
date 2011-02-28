SUBROUTINE matr_mesh(matr,nmat,VB,nvb,IT,nit,dif)

IMPLICIT NONE
	integer ::nit, IT(3,nit),ittmp(3),i,nvb,nmat
	real*8 :: point(1,1:3),clost(1,1:3),VB(3,nvb),normals(3,3),dist(1,3),dif(nmat),DAT(nit,12),matr(nmat,3),diff
 
 call updateSearch(VB,nvb,IT,nit,DAT)
 
 do i = 1,nmat

    point(1,1:3) = matr(i,:)
    
    call pt_upmesh(point,DAT,nit,clost,diff)
  ! dist(1,1:3) = (clostmp(1,1:3)-point(1,1:3))
   ! dif = dot_product(dist(1,1:3),dist(1,1:3))
    matr(i,:) = clost(1,1:3)
    dif(i)=diff
 end do

END SUBROUTINE matr_mesh
