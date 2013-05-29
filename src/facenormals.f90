SUBROUTINE facenormals(VB,nvb,IT,nit,O,normals)
  
  IMPLICIT NONE
  integer :: nvb,i,j,t,nit,IT(3,nit),O
  real*8 :: VB(3,nvb),l,ntmp(3),tmp0(1,1:3),tmp1(1,1:3),normals(3,nit)
  real*8 :: tmp2(1,1:3)
  
 
 ! pi = ACOS(0.0)
  !PI=4.D0*DATAN(1.D0)

  data ntmp / 0,0,0 /
  do i = 1,nit
     tmp0(1,1:3) = VB(1:3,IT(2,i))-VB(1:3,IT(1,i))
     tmp1(1,1:3) = VB(1:3,IT(3,i))-VB(1:3,IT(1,i))
     call crossp(tmp0,tmp1,ntmp)
     normals(1:3,i) = ntmp
  end do
  do i = 1,nit
     call normalize(normals(1:3,i))
  end do
contains
  subroutine veclen (x,l)
    real*8 :: x(3),l
    l=sqrt(sum(x**2))
    
  end subroutine veclen
  subroutine normalize (x)
    real*8 :: x(3),l
    call veclen(x,l)	
    if (l > 0) then
    x = x/l
    end if 
  end subroutine normalize
END SUBROUTINE facenormals


