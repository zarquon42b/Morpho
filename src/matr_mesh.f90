SUBROUTINE matr_mesh(matr,nmat,VB,nvb,IT,nit,dif)

IMPLICIT NONE
	integer ::nit, IT(3,nit),ittmp(3),i,nvb,nmat
	real :: point(3),clost(3),VB(3,nvb),normals(3,3),dif(nmat),DAT(nit,12),matr(nmat,3),diff
 
 call updateSearch(VB,nvb,IT,nit,DAT)
 
 do i = 1,nmat

    point(:) = matr(i,:)
    
    call pt_upmesh(point,DAT,nit,clost,diff)
    matr(i,:) = clost(:)
    dif(i)=diff
 end do

END SUBROUTINE matr_mesh
