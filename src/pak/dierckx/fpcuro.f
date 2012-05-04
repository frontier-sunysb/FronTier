      subroutine fpcuro(a,b,c,d,x,n)
      implicit none
c  subroutine fpcuro finds the double precision zeros of a cubic polynomial
c  p(x) = a*x**3+b*x**2+c*x+d.
c
c  calling sequence:
c     call fpcuro(a,b,c,d,x,n)
c
c  input parameters:
c    a,b,c,d: double precision values, containing the coefficients of p(x).
c
c  output parameters:
c    x      : double precision array,length 3, which contains the double precision zeros of p(x)
c    n      : integer, giving the number of double precision zeros of p(x).
c  ..
c  ..scalar arguments..
      double precision a,b,c,d
      integer n
c  ..array argument..
      double precision x(3)
c  ..local scalars..
      integer i
      double precision a1,b1,c1,df,disc,d1,e3,f,four,half,ovfl,pi3,p3,
     * q,r,step,tent,three,two,u,u1,u2,y
c  ..function references..
      double precision abs,dmax1,atan,datan2,cos,sign,dsqrt
c  set constants
      two = 0.2d+01
      three = 0.3d+01
      four = 0.4d+01
      ovfl =0.1d+05
      half = 0.5d+0
      tent = 0.1d+0
      e3 = tent/0.3d0
      pi3 = atan(0.1d+01)/0.75d0
      a1 = abs(a)
      b1 = abs(b)
      c1 = abs(c)
      d1 = abs(d)
c  test whether p(x) is a third degree polynomial.
      if(dmax1(b1,c1,d1).lt.a1*ovfl) go to 300
c  test whether p(x) is a second degree polynomial.
      if(dmax1(c1,d1).lt.b1*ovfl) go to 200
c  test whether p(x) is a first degree polynomial.
      if(d1.lt.c1*ovfl) go to 100
c  p(x) is a constant function.
      n = 0
      go to 800
c  p(x) is a first degree polynomial.
 100  n = 1
      x(1) = -d/c
      go to 500
c  p(x) is a second degree polynomial.
 200  disc = c*c-four*b*d
      n = 0
      if(disc.lt.0.) go to 800
      n = 2
      u = dsqrt(disc)
      b1 = b+b
      x(1) = (-c+u)/b1
      x(2) = (-c-u)/b1
      go to 500
c  p(x) is a third degree polynomial.
 300  b1 = b/a*e3
      c1 = c/a
      d1 = d/a
      q = c1*e3-b1*b1
      r = b1*b1*b1+(d1-b1*c1)*half
      disc = q*q*q+r*r
      if(disc.gt.0.) go to 400
      u = dsqrt(abs(q))
      if(r.lt.0.) u = -u
      p3 = datan2(dsqrt(-disc),abs(r))*e3
      u2 = u+u
      n = 3
      x(1) = -u2*cos(p3)-b1
      x(2) = u2*cos(pi3-p3)-b1
      x(3) = u2*cos(pi3+p3)-b1
      go to 500
 400  u = dsqrt(disc)
      u1 = -r+u
      u2 = -r-u
      n = 1
      x(1) = sign(abs(u1)**e3,u1)+sign(abs(u2)**e3,u2)-b1
c  apply a newton iteration to improve the accuracy of the roots.
 500  do 700 i=1,n
        y = x(i)
        f = ((a*y+b)*y+c)*y+d
        df = (three*a*y+two*b)*y+c
        step = 0.
        if(abs(f).lt.abs(df)*tent) step = f/df
        x(i) = y-step
 700  continue
 800  return
      end
