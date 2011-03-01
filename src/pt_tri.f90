SUBROUTINE pt_tri(point,VBvec,clost)

IMPLICIT NONE

 real*8 :: point(3),clost(3)
 real*8 :: d,e,f,det,s,t,numer,denom,dv(3)
 real*8 :: VBvec(12)
 real*8 :: B1(3),e0(3),e1(3),a,b,c

 B1 = VBvec(1:3)
 dv = point(1:3)-B1(1:3)
 e0 = VBvec(4:6)
 e1 = VBvec(7:9)
	!data ntmp(1,1:3) / 0,0,0 /
 a = VBvec(10)
 b = VBvec(11)
 c = VBvec(12)
 d = dot_product(e0(:),dv(:))
 e = dot_product(e1(:),dv(:))
 f = dot_product(dv(:),dv(:))

 det = a*c - b*b
 s = b*2
 t=b*d - a*e
 !if (s+t <= det) then
 !   if (s < 0) then
  !     if (t < 0) then
  !        !region 4 

 det = 1/det
 s = s*det
 t= t*det
 numer = c+d-b-d
 if (numer <=0) then
    s = 0
 else
    denom = a -2 * b +c
    if (numer >= denom) then
       s = numer
    else
    s = (numer/denom)
 end if
 end if
 t = 1-s
 

 clost(:) = B1(:)+ s*e0(:) + t*e1(:)
 

END SUBROUTINE pt_tri
