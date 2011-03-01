SUBROUTINE pt_angmesh(point,normal,DAT,ndat,clost,dif,fptr,rhotol,tarnorm)

IMPLICIT NONE
integer ::ndat,i,fptr,ittmp(3)
real :: point(3),clost(3),clostmp(3),dist(3),dif,dif_old,DAT(ndat,15),clostang(3)
real :: vbtmp(12),normal(3),rhotol,rho,tarnorm(3,ndat),lnorm,normtmp(3)

 clostang(:) = (/9999,9999,9999/)
 dif_old=1e10
 
 do i = 1,ndat
   
    vbtmp(:) = DAT(i,1:12)
    ittmp(:) = int(DAT(i,13:15))
    call pt_tri(point,vbtmp,clostmp)
    dist(:) = (clostmp(:)-point(:))
    dif = dot_product(dist(:),dist(:))
    
    if (dif <= dif_old) then

       clost = clostmp

       !calculate normal from surrounding vertices
       normtmp(:) = tarnorm(:,ittmp(1))+tarnorm(:,ittmp(2))+tarnorm(:,ittmp(3))
       lnorm = sqrt(dot_product(normtmp(:),normtmp(:)))
       normtmp = normtmp/lnorm
      ! call angcal(normtmp,3,normal,3,rho)
      ! if (rho <= rhotol) then
      !    dif_old = dif
      !    clostang = clostmp
     !     fptr = i
     ! end if
   end if
      
 end do
 if  (.NOT.(clostang(1) .eq. 9999)) then
    clost = clostang
 end if
END SUBROUTINE pt_angmesh

