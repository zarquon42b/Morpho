SUBROUTINE pt_upmesh(point,DAT,ndat,clost,diff,fptr,region,met)

IMPLICIT NONE
integer :: ndat,i,fptr,region,regtmp,method
real*8 :: point(3),clost(3),clostmp(3),dist(3),dif,tmpdist
real*8 :: dif_old,DAT(ndat,13), vbtmp(13),diff,sqdist, normal(3),difvec(3)
real*8 :: checkclost(3)
integer :: met
clostmp(:) = (/9999,9999,9999/)
dif_old=1e10
!method=0
!if (present(met)) then
!   method=met
!methodend if
!write (*,*),met
do i = 1,ndat
   
   vbtmp(1:13) = DAT(i,1:13)
   
   !! check if distance to triangle plane exceeds dif_old
   !!call  pt_triplane(point,vbtmp,checkclost,tmpdist)
   !!if (tmpdist < dif_old) then
      call pt_tri(point,vbtmp,clostmp,regtmp,sqdist)
      if (met == 1 .and. regtmp /= 0 ) then
         call pt_triplane(point,vbtmp(1:12),checkclost,tmpdist)
         sqdist = (sqrt(tmpdist) + sqrt(sum((checkclost-clostmp)**2)))**2
      end if
         
      !dist(:) = sum(clostmp(:)-point(:))
      !dif = dot_product(dist(:),dist(:))
      
      if (sqdist <= dif_old) then
         diff = abs(sqdist)
         diff = sqrt(diff)
         dif_old = sqdist
         clost = clostmp
         fptr = i
         region = regtmp
      end if
  !! end if
end do



END SUBROUTINE pt_upmesh
