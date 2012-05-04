      double precision function evapol(tu,nu,tv,nv,c,rad,x,y)
      implicit none
c  function program evacir evaluates the function f(x,y) = s(u,v),
c  defined through the transformation
c      x = u*rad(v)*cos(v)    y = u*rad(v)*sin(v)
c  and where s(u,v) is a bicubic spline ( 0<=u<=1 , -pi<=v<=pi ), given
c  in its standard b-spline representation.
c
c  calling sequence:
c     f = evapol(tu,nu,tv,nv,c,rad,x,y)
c
c  input parameters:
c   tu    : double precision array, length nu, which contains the position of the
c           knots in the u-direction.
c   nu    : integer, giving the total number of knots in the u-direction
c   tv    : double precision array, length nv, which contains the position of the
c           knots in the v-direction.
c   nv    : integer, giving the total number of knots in the v-direction
c   c     : double precision array, length (nu-4)*(nv-4), which contains the
c           b-spline coefficients.
c   rad   : double precision function subprogram, defining the boundary of the
c           approximation domain. must be declared external in the
c           calling (sub)-program
c   x,y   : double precision values.
c           before entry x and y must be set to the co-ordinates of
c           the point where f(x,y) must be evaluated.
c
c  output parameter:
c   f     : real
c           on exit f contains the value of f(x,y)
c
c  other subroutines required:
c    bispev,fpbisp,fpbspl
c
c  references :
c    de boor c : on calculating with b-splines, j. approximation theory
c                6 (1972) 50-62.
c    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
c                applics 10 (1972) 134-149.
c    dierckx p. : curve and surface fitting with splines, monographs on
c                 numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1989
c
c  ..scalar arguments..
      integer nu,nv
      double precision x,y
c  ..array arguments..
      double precision tu(nu),tv(nv),c((nu-4)*(nv-4))
c  ..user specified function
      double precision rad
c  ..local scalars..
      integer ier
      double precision u,v,r,f,one,dist
c  ..local arrays
      double precision wrk(8)
      integer iwrk(2)
c  ..function references
      double precision datan2,dsqrt
c  ..
c  calculate the (u,v)-coordinates of the given point.
      one = 1
      u = 0.
      v = 0.
      dist = x**2+y**2
      if(dist.le.0.) go to 10
      v = datan2(y,x)
      r = rad(v)
      if(r.le.0.) go to 10
      u = dsqrt(dist)/r
      if(u.gt.one) u = one
c  evaluate s(u,v)
  10  call bispev(tu,nu,tv,nv,c,3,3,u,1,v,1,f,wrk,8,iwrk,2,ier)
      evapol = f
      return
      end

