SUBROUTINE pt_angmesh(point,normal,DAT,ndat,itnorm,clost,dif,fptr,rhotol)

  IMPLICIT NONE
  integer :: ndat,i,fptr,region,regtmp,ittmp(3)
  real*8 :: point(3),clost(3),clostmp(3),clostang(3),dist(3),dif
  real*8 :: normal(3),lnorm,normtmp(3),itnorm(ndat,3)
  real*8 :: dif_old,DAT(ndat,12), vbtmp(12),diff,sqdist
  real*8 :: rho,rhotol
  !real*8, target :: DAT(ndat,12)
  !real*8, pointer :: vbtmp(:)

  !allocate(vbtmp(12))
  
  dif_old=1e10
  clostang(:) = (/9999,9999,9999/)
  dif_old=1e10
  
  do i = 1,ndat

     vbtmp(1:12) = DAT(i,1:12)

     call pt_tri(point,vbtmp,clostmp,regtmp,sqdist)


     if (sqdist <= dif_old) then

        clost = clostmp
        !calculate normal from surrounding vertices

        !call angcal(itnorm(i,:),3,normal,3,rho)
       ! if (rho <= rhotol) then
        !   diff = sqdist
       !    dif_old = sqdist
       !    clostang = clostmp
           fptr = i
          region = regtmp
      ! end if
       end if

    end do

   ! if  (.NOT.(clostang(1) .eq. 9999)) then
   !  clost = clostang
 ! end if


END SUBROUTINE pt_angmesh

