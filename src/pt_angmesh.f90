SUBROUTINE pt_angmesh(point,normal,DAT,ndat,itnorm,clost,diff,fptr,rhotol,toldist)

  IMPLICIT NONE
  integer :: ndat,i,fptr,fptr2,region,regtmp,ittmp(3)
  real*8 :: point(3),clost(3),clostmp(3),clostang(3),dist(3),dif
  real*8 :: normal(3),lnorm,normtmp(3),itnorm(ndat,3),itnormtmp(3)
  real*8 :: dif_old,DAT(ndat,12), vbtmp(12),diff,diff2,sqdist,toldist
  real*8 :: rho,rhotol,dif_old2
  !real*8, target :: DAT(ndat,12)
  !real*8, pointer :: vbtmp(:)

  !allocate(vbtmp(12))
  !rhotol = 0.6
  dif_old = 1e10
  dif_old2 = 1e10
  clostang(:) = (/9999,9999,9999/)
  !print *, itnorm(1,:)
  
  do i = 1,ndat
     itnormtmp(:) = itnorm(i,:)
     vbtmp(1:12) = DAT(i,1:12)

     call pt_tri(point,vbtmp,clostmp,regtmp,sqdist)


     if (sqdist <= dif_old2) then
      
        clost = clostmp
        diff2 = sqdist
        dif_old2 = sqdist
        fptr2 = i
        
       ! calculate normal from surrounding vertices
        
     end if
     
     if (sqdist <= dif_old ) then
        call angcal(itnormtmp,3,normal,3,rho)
        if (rho <= rhotol .AND. sqdist <= toldist) then
           diff = sqdist
           dif_old = sqdist
           clostang = clostmp
           fptr = i
          ! region = regtmp
           ! end if
        end if
     
     end if
  end do

  if  (.NOT.(clostang(1) .eq. 9999)) then
        clost = clostang
     else
        diff = diff2
        fptr = fptr2
     end if


   END SUBROUTINE pt_angmesh

