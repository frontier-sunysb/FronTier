      subroutine drotm (n,dx,incx,dy,incy,dparam)
c
c     apply the modified givens transformation, h, to the 2 by n matrix
c
c     (dx**t) , where **t indicates transpose. the elements of dx are in
c     (dy**t)
c
c     dx(lx+i*incx), i = 0 to n-1, where lx = 1 if incx .ge. 0, else
c     lx = (-incx)*n, and similarly for sy using ly and incy.
c     with dparam(1)=dflag, h has one of the following forms..
c
c     dflag=-1.d0     dflag=0.d0        dflag=1.d0     dflag=-2.d0
c
c       (dh11  dh12)    (1.d0  dh12)    (dh11  1.d0)    (1.d0  0.d0)
c     h=(          )    (          )    (          )    (          )
c       (dh21  dh22),   (dh21  1.d0),   (-1.d0 dh22),   (0.d0  1.d0).
c     see drotmg for a description of data storage in dparam.
c
      double precision dflag,dh12,dh22,dx,two,z,dh11,dh21,
     1 dparam,dy,w,zero
      dimension dx(1),dy(1),dparam(5)
      data zero,two/0.d0,2.d0/
c
      dflag=dparam(1)
      if(n .le. 0 .or.(dflag+two.eq.zero)) go to 140
          if(.not.(incx.eq.incy.and. incx .gt.0)) go to 70
c
               nsteps=n*incx
               if(dflag) 50,10,30
   10          continue
               dh12=dparam(4)
               dh21=dparam(3)
                    do 20 i=1,nsteps,incx
                    w=dx(i)
                    z=dy(i)
                    dx(i)=w+z*dh12
                    dy(i)=w*dh21+z
   20               continue
               go to 140
   30          continue
               dh11=dparam(2)
               dh22=dparam(5)
                    do 40 i=1,nsteps,incx
                    w=dx(i)
                    z=dy(i)
                    dx(i)=w*dh11+z
                    dy(i)=-w+dh22*z
   40               continue
               go to 140
   50          continue
               dh11=dparam(2)
               dh12=dparam(4)
               dh21=dparam(3)
               dh22=dparam(5)
                    do 60 i=1,nsteps,incx
                    w=dx(i)
                    z=dy(i)
                    dx(i)=w*dh11+z*dh12
                    dy(i)=w*dh21+z*dh22
   60               continue
               go to 140
   70     continue
          kx=1
          ky=1
          if(incx .lt. 0) kx=1+(1-n)*incx
          if(incy .lt. 0) ky=1+(1-n)*incy
c
          if(dflag)120,80,100
   80     continue
          dh12=dparam(4)
          dh21=dparam(3)
               do 90 i=1,n
               w=dx(kx)
               z=dy(ky)
               dx(kx)=w+z*dh12
               dy(ky)=w*dh21+z
               kx=kx+incx
               ky=ky+incy
   90          continue
          go to 140
  100     continue
          dh11=dparam(2)
          dh22=dparam(5)
               do 110 i=1,n
               w=dx(kx)
               z=dy(ky)
               dx(kx)=w*dh11+z
               dy(ky)=-w+dh22*z
               kx=kx+incx
               ky=ky+incy
  110          continue
          go to 140
  120     continue
          dh11=dparam(2)
          dh12=dparam(4)
          dh21=dparam(3)
          dh22=dparam(5)
               do 130 i=1,n
               w=dx(kx)
               z=dy(ky)
               dx(kx)=w*dh11+z*dh12
               dy(ky)=w*dh21+z*dh22
               kx=kx+incx
               ky=ky+incy
  130          continue
  140     continue
          return
          end
