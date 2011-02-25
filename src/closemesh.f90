SUBROUTINE closemesh(point,VB,nvb,IT,nit,normals,clost,dif)

IMPLICIT NONE
	integer ::nit, IT(nit,3),ittmp(1:3),i,nvb,j
	real*8 :: point(1,1:3),clost(1,1:3),clostmp(1,1:3),VB(nvb,3),normals(3,3),vbtmp(3,3),dist(1,1:3),dif,dif_old
 ittmp =  (/1,2,3/)
 clostmp(1,1:3) = (/0,0,0/)
 dif_old = 1000000
 do i = 1,2
    ittmp(1:3) = IT(1,1:3)
    vbtmp (1:3,1:3) = VB(ittmp(1:3),1:3)
  
    !vbtmp(1:3,1:3) = VB(ittmp,1:3)
    
    !call closest(point,vbtmp,normals,clostmp)

    !dist(1,1:3) = (clost(1,1:3)-point(1,1:3))
    !dif = dot_product(dist(1,1:3),dist(1,1:3))
    !if (dif < dif_old) then
     !  dif_old = dif
     !  clost = clostmp
    !end if
       

 
 end do

END SUBROUTINE closemesh
