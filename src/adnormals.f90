SUBROUTINE adnormals(VB,nvb,IT,nit,O,normals,angweight)
  
  IMPLICIT NONE
  integer :: nvb,i,j,t,nit,IT(3,nit),O
  real*8 :: VB(3,nvb),l,ntmp(1,1:3),tmp0(1,1:3),tmp1(1,1:3),normals(4,nvb),co
  real*8 ::angtmp(3),pi,tmp2(1,1:3)
  logical :: angweight

 ! pi = ACOS(0.0)
  !PI=4.D0*DATAN(1.D0)

  data ntmp(1,1:3) / 0,0,0 /
  do i = 1,nit
     tmp0(1,1:3) = VB(1:3,IT(2,i))-VB(1:3,IT(1,i))
     tmp1(1,1:3) = VB(1:3,IT(3,i))-VB(1:3,IT(1,i))
   if (angweight .eqv. .true.) then 
     tmp2(1,1:3) = VB(1:3,IT(2,i))-VB(1:3,IT(3,i))
     call angcal (tmp0,3,tmp1,3,angtmp(1))
     call angcal (-tmp1,3,tmp2,3,angtmp(2))
     angtmp(3) = 3.141592653589793239  -angtmp(1)-angtmp(2)
  end if

     call crossp(tmp0,tmp1,ntmp)
     !call veclen(ntmp,l)
     !if (l==0) then	
     !   ntmp(1,1:3) = (/0 , 0 , 0/)
        
     !else
     !   ntmp = ntmp/l
	
     !end if
		
     do j = 1,3
        co =sum(normals(1:3,IT(j,i))*ntmp(1,1:3))
        if (co<0) then
           if (angweight .eqv. .false.) then 
              normals(1:3,IT(j,i)) = normals(1:3,IT(j,i))-ntmp(1,1:3)
              !call normalize(normals(1:3,IT(j,i)),l)
              !normals(4,IT(j,i)) = normals(4,IT(j,i)) + 1
           else
              normals(1:3,IT(j,i)) = normals(1:3,IT(j,i))-ntmp(1,1:3)*angtmp(j)
           end if
        else 
           if (angweight .eqv. .false.) then 
              normals(1:3,IT(j,i)) = normals(1:3,IT(j,i))+ntmp(1,1:3)
              !normals(1:3,IT(j,i)) = normals(1:3,IT(j,i))
              !call normalize(normals(1:3,IT(j,i)),l)
              !normals(4,IT(j,i)) = normals(4,IT(j,i)) + 1
           else
              normals(1:3,IT(j,i)) = normals(1:3,IT(j,i))+ntmp(1,1:3)*angtmp(j)
           end if
        end if
     end do
  end do
  do i = 1,nvb
    !  if (angweight .eqv. .false.) then 
    ! normals(1:4,i) = normals(1:4,i)/normals(4,i)
!end if
     !call veclen(normals,l)
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
END SUBROUTINE adnormals


