SUBROUTINE rho_mesh(matr,nmat,normals,VB,nvb,IT,nit,tarnorm,dif,fptr,rhotol)

IMPLICIT NONE
integer :: nit, IT(3,nit),ittmp(3),i,nvb,nmat,fptr
real :: point(3),clost(3),VB(3,nvb),dif(nmat),DAT(nit,15),matr(nmat,3)
real :: diff,rhotol,normals(3,nvb),normtmp(3),tarnorm(3,nvb)
 
 call updateSearch(VB,nvb,IT,nit,DAT)
 
 do i = 1,nmat
    
    point(:) = matr(i,:)
    normtmp(:) = normals(:,i)
    call pt_angmesh(point,DAT,nit,clost,diff,fptr,rhotol,tarnorm)
    matr(i,:) = clost(:)
    dif(i)=diff
 end do

END SUBROUTINE rho_mesh
