SUBROUTINE minClost(X,x1,Y,y1,out)
  
  IMPLICIT NONE
  integer :: 	x1,y1,i,j,out(x1)
  real*8 :: 	X(x1,3), Y(y1,3),tmp, tmpvec(3), dist
  
 do i = 1,x1
    dist=1e10
    do j = 1,y1
       tmpvec = X(i,1:3) - Y(j,1:3)
       tmp = DOT_PRODUCT(tmpvec,tmpvec)
       if (tmp < dist) then
          dist = tmp
          out(i) = j
       end if
    end do
 end do
 
END SUBROUTINE minClost
