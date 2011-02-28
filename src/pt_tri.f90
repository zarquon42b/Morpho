SUBROUTINE pt_tri(point,VBvec,clost)

IMPLICIT NONE
	!integer ::
	real :: point(3),clost(3)
 real ::d,e,f,det,s,t,numer,denom,dv(3)
 real :: VBvec(12)
 real :: B1(3),e0(3),e1(3),a,b,c

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
 det = 1/det
 s = s*det
 t= t*det
 numer = c+d-b-d
 if (numer <=0) then
    s = 0
 else
    denom = a -2 * b +c
    s = (numer/denom)
    
 end if
 t = 1-s
 

 clost(:) = B1(:)+ s*e0(:) + t*e1(:)
 
 !c = matmul(transpose(e1),e1)
 !do i = 1,N

END SUBROUTINE pt_tri
