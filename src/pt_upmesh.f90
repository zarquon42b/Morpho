SUBROUTINE pt_upmesh(point,DAT,ndat,clost,dif)

IMPLICIT NONE
	integer ::ndat,i
	real*8 :: point(1,1:3),clost(1,1:3),clostmp(1,1:3),dist(1,3),dif,dif_old,DAT(ndat,12),vbtmp(12)
 clostmp(1,1:3) = (/9999,9999,9999/)
 dif_old=1e10
  
 do i = 1,ndat
   
    vbtmp = DAT(i,1:12)
    
    call pt_tri(point,vbtmp,clostmp)
    dist(1,1:3) = (clostmp(1,1:3)-point(1,1:3))
    dif = dot_product(dist(1,1:3),dist(1,1:3))
   
    if (dif <= dif_old) then
      dif_old = dif
      clost = clostmp
      
   end if
      

 
 end do

END SUBROUTINE pt_upmesh
