SUBROUTINE updateSearch(VB,nvb,IT,nit,DAT)
!! sort data to use for closest point on triangle search

IMPLICIT NONE
	integer :: nit,IT(3,nit),ittmp(3),i,nvb
	real*8 :: VB(3,nvb),DAT(nit,13),vbtmp(3,3),e0(3),e1(3),B(3)
! real*8, pointer :: e0(:),e1(:),B(:)
 
 do i=1,nit
    
    ittmp(:) = IT(:,i)
    vbtmp (:,:) = transpose(VB(:,ittmp))
    !print *,vbtmp
    
    DAT(i,1:3) = vbtmp(1,1:3)!B1
    B(:) = DAT(i,1:3)
    DAT(i,4:6) = vbtmp(2,1:3)-vbtmp(1,1:3)!e0
    e0(1:3) = DAT(i,4:6)
    DAT(i,7:9) = vbtmp(3,1:3)-vbtmp(1,1:3)!e1
    e1(1:3) = DAT(i,7:9)
    DAT(i,10) = dot_product(e0(:),e0(:))!a
    DAT(i,11) = dot_product(e0(:),e1(:))!b
    DAT(i,12) = dot_product(e1(:),e1(:))!c
    DAT(i,13) = abs(DAT(i,10)*DAT(i,12) - DAT(i,11)**2)

end do
 

END SUBROUTINE updateSearch
