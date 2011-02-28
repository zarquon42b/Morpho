SUBROUTINE pt_mesh(point,VB,nvb,IT,nit,clost,dif)

IMPLICIT NONE
	integer ::nit, IT(3,nit),ittmp(3),i,nvb,j
	real*8 :: point(1,1:3),clost(1,1:3),clostmp(1,1:3),VB(3,nvb),normals(3,3),dist(1,3),dif,dif_old,DAT(nit,12),vbtmp(12)
 clostmp(1,1:3) = (/9999,9999,9999/)
 dif_old=1e10
 call updateSearch(VB,nvb,IT,nit,DAT)
 
 do i = 1,nit
   
    vbtmp = DAT(i,1:12)
    !print *,vbtmp  
    call pt_tri(point,vbtmp,clostmp)
    dist(1,1:3) = (clostmp(1,1:3)-point(1,1:3))
    dif = dot_product(dist(1,1:3),dist(1,1:3))
   !print *, dif
    !print * ,ittmp
    !print *, clostmp
    if (dif <= dif_old) then
      dif_old = dif
      clost = clostmp
      
   end if
      

 
 end do

END SUBROUTINE pt_mesh
