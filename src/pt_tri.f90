SUBROUTINE pt_tri(point,VBvec,clost,region,sqdist)
!!! calculate distance between point and triangle
IMPLICIT NONE
integer :: region
 real*8 :: point(3),clost(3)
 real*8 :: det,s,t,numer,denom,dv(3),invDet,tmp0,tmp1
 real*8 :: VBvec(1:13)
 real*8 :: B(3),e0(3),e1(3)
 real*8 :: a00,a01,a11,b0,b1,c,sqdist

 B(:) = VBvec(1:3)
 dv(1:3) = B(1:3) - point(1:3)
 e0(1:3) = VBvec(4:6)
 e1(1:3) = VBvec(7:9)
 a00 = VBvec(10)
 a01 = VBvec(11)
 a11 = VBvec(12)
 b0 = dot_product(e0(1:3),dv(1:3))
 b1 = dot_product(e1(1:3),dv(1:3))
 c = dot_product(dv(1:3),dv(1:3))
 !det = abs(a00*a11 - a01*a01)
 det =  VBvec(13)
 s = a01*b1- a11*b0
 t = a01*b0 - a00*b1

 if (s+t <= det) then
    
    if (s < 0) then
       if (t < 0) then

          !region 4 begin
          region = 4

          if (b0 < 0) then
             t = 0
             if ( -b0 >= a00) then
                s = 1
                sqdist = a00 + 2*b0 + c
             else
                s = -b0/a00
                sqdist = b0*s + c
             end if
         
          else
             s = 0
             if ( b1 >= 0 ) then
                t=0
                sqdist = c
             else if (-b1 >= a11) then
                t = 1
                sqdist = a11 + 2*b1 + c
             else
                t = -b1/a11
                sqdist = b1*t +c
             end if
          end if
       else
          !region 3 begin
          region = 3
          s = 0
          if (b1 >=0) then
             t=0
             sqdist = c
          else if (-b1 >= a11) then
             t = 1
             sqdist = a11+2*b1+c
          else
             t = -b1/a11
             sqdist = b1*t+c
          end if
          !region3 end
       end if
       

    else if (t < 0) then
       !region 5 begin
       region = 5
       t = 0
       if (b0 >= 0) then
          s = 0
          sqdist = c
       else if (-b0 >=a00) then
          s = 1
          sqdist = a00 + 2*b0 + c
       else 
          s = -b0/a00
          sqdist = b0*s + c

       end if
       !region 5 end

    else        
       !region 0 begin
       region = 0
       invDet = 1/det
       s = s*invDet
       t=  t*invDet
       sqdist = s*(a00*s + a01*t + 2*b0) + t*(a01*s + a11*t + 2*b1) + c
       ! region 0 end
    end if


 else
    if (s < 0) then
       !region 2 begin
       region = 2
       tmp0 = a01 + b0
       tmp1 = a11 + b1
       if (tmp1 > tmp0) then
          numer = tmp1 - tmp0
          denom = a00 - 2*a01 + a11

          if (numer >= denom) then
             s = 1
             t = 0
             sqdist = a00 + 2*b0 + c

          else               
             s = numer/denom
             t = 1 - s
             sqdist = s*(a00*s + a01*t + 2*b0) + t*(a01*s + a11*t +2*b1) + c
          end if

       else 
          s = 0
          if (tmp1 <= 0) then
             t = 1
             sqdist = a11+2*b1+c
          else if (b1 >= 0) then
             t = 0
             sqdist = c
          else
             t = -b1/a11
             sqdist = b1*t + c
          end if
          !region 2 end
       end if
    else if (t < 0) then
       ! region 6 begin
       region = 6
       tmp0 = a01 + b1
       tmp1 = a00 + b0
       if (tmp1 > tmp0) then
          numer = tmp1 - tmp0
          denom = a00 - 2*a01 + a11
          if (numer >= denom) then                             
             t = 1
             s = 0
             sqdist = a11 + 2*b1 + c
          else 
             t = numer/denom
             s = 1 - t
             sqdist = s*(a00*s + a01*t + 2*b0) + t*(a01*s + a11*t + 2*b1) + c
          end if

       else
          t = 0
          if (tmp1 <= 0) then
             s = 1
             sqdist = a00 + 2*b0 + c
          else if ( b0 >=0) then
             s = 0
             sqdist = c
          else
             s = -b0/a00
             sqdist = b0*s +c
          end if
       end if
       else
          numer = a11 + b1 - a01 - b0
          if (numer <= 0) then
             s = 0
             t = 1 
             sqdist = a11+ 2*b1 + c
          else
             denom = a00 - 2*a01 + a11;
             if (numer >= denom) then

             s =1
             t = 0
             sqdist = a00 + 2*b0 + c

          else

             s = numer/denom
             t = 1 - s
             sqdist = s*(a00*s + a01*t + 2*b0) + t*(a01*s + a11*t + 2*b1) + c

          end if
          
       end if
    end if





 end if!// closes first if


 
 !sqdist = sqdist+1
 clost = B + s*e0 + t*e1
END SUBROUTINE pt_tri


subroutine pt_triplane(point,VBvec,clost,sqdist)
!!! calculate distance between point and triangle plane
!!! updates clost and sqdist (closest point on plane and squared distance)
IMPLICIT NONE
integer :: region
real*8 :: point(3),clost(3)
real*8 :: det,s,t,numer,denom,dv(3),invDet,tmp0,tmp1,normal(3)
real*8 :: VBvec(1:12),difvec(3)
real*8 :: B(3),e0(3),e1(3),p0p0(3)
real*8 :: a00,a01,a11,b0,b1,c,sqdist,p0p0dist,alpha

B(:) = VBvec(1:3)
dv(1:3) = B(1:3) - point(1:3)
e0(1:3) = VBvec(4:6)
e1(1:3) = VBvec(7:9)
call crossp(e0,e1,normal)
normal = normal/sqrt(dot_product(normal,normal))
difvec = point-B
difvec = difvec/sqrt(dot_product(difvec,difvec))
alpha = dot_product(difvec,normal)
p0p0dist = sqrt(dot_product(B -point,B -point))*alpha
p0p0 = -p0p0dist*normal
clost = point + p0p0
sqdist = p0p0dist**2
END SUBROUTINE pt_triplane


