SUBROUTINE pt_upmesh(point,DAT,ndat,clost,diff,fptr,region)

IMPLICIT NONE
integer :: ndat,i,fptr,region,regtmp
real*8 :: point(3),clost(3),clostmp(3),dist(3),dif
real*8 :: dif_old,DAT(ndat,12), vbtmp(12),diff,sqdist, normal(3),difvec(3)
 !real*8, target :: DAT(ndat,12)
 !real*8, pointer :: vbtmp(:)

 !allocate(vbtmp(12))
 clostmp(:) = (/9999,9999,9999/)
 dif_old=1e10
 
 do i = 1,ndat
   
    vbtmp(1:12) = DAT(i,1:12)

!! check if distance to triangle plane exceeds dif_old
!!    call crossp(DAT(i,4:6),DAT(i,7:9),normal)
!!    normal = normal/sqrt(DOT_PRODUCT(normal,normal))
!!    difvec = DAT(i,1:3) - point 
!!    sqdist = sqrt(dot_product(normal,difvec))
    
!!    if (abs(sqdist) < dif_old) then
       call pt_tri(point,vbtmp,clostmp,regtmp,sqdist)
       
       !dist(:) = sum(clostmp(:)-point(:))
       !dif = dot_product(dist(:),dist(:))
       
       if (sqdist <= dif_old) then
          diff = sqrt(sqdist)
          dif_old = sqdist
          clost = clostmp
          fptr = i
          region = regtmp
          !print *,sqdist
       end if
!!    end if
 end do
 
END SUBROUTINE pt_upmesh
