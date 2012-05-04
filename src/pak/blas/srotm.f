      subroutine srotm (n,sx,incx,sy,incy,sparam)
c
c     apply the modified givens transformation, h, to the 2 by n matrix
c
c     (sx**t) , where **t indicates transpose. the elements of sx are in
c     (dx**t)
c
c     sx(lx+i*incx), i = 0 to n-1, where lx = 1 if incx .ge. 0, else
c     lx = (-incx)*n, and similarly for sy using using ly and incy.
c     with sparam(1)=sflag, h has one of the following forms..
c
c     sflag=-1.e0     sflag=0.e0        sflag=1.e0     sflag=-2.e0
c
c       (sh11  sh12)    (1.e0  sh12)    (sh11  1.e0)    (1.e0  0.e0)
c     h=(          )    (          )    (          )    (          )
c       (sh21  sh22),   (sh21  1.e0),   (-1.e0 sh22),   (0.e0  1.e0).
c     see  srotmg for a description of data storage in sparam.
c
      dimension sx(1),sy(1),sparam(5)
      data zero,two/0.e0,2.e0/
c
      sflag=sparam(1)
      if(n .le. 0 .or.(sflag+two.eq.zero)) go to 140
          if(.not.(incx.eq.incy.and. incx .gt.0)) go to 70
c
               nsteps=n*incx
               if(sflag) 50,10,30
   10          continue
               sh12=sparam(4)
               sh21=sparam(3)
                    do 20 i=1,nsteps,incx
                    w=sx(i)
                    z=sy(i)
                    sx(i)=w+z*sh12
                    sy(i)=w*sh21+z
   20               continue
               go to 140
   30          continue
               sh11=sparam(2)
               sh22=sparam(5)
                    do 40 i=1,nsteps,incx
                    w=sx(i)
                    z=sy(i)
                    sx(i)=w*sh11+z
                    sy(i)=-w+sh22*z
   40               continue
               go to 140
   50          continue
               sh11=sparam(2)
               sh12=sparam(4)
               sh21=sparam(3)
               sh22=sparam(5)
                    do 60 i=1,nsteps,incx
                    w=sx(i)
                    z=sy(i)
                    sx(i)=w*sh11+z*sh12
                    sy(i)=w*sh21+z*sh22
   60               continue
               go to 140
   70     continue
          kx=1
          ky=1
          if(incx .lt. 0) kx=1+(1-n)*incx
          if(incy .lt. 0) ky=1+(1-n)*incy
c
          if(sflag)120,80,100
   80     continue
          sh12=sparam(4)
          sh21=sparam(3)
               do 90 i=1,n
               w=sx(kx)
               z=sy(ky)
               sx(kx)=w+z*sh12
               sy(ky)=w*sh21+z
               kx=kx+incx
               ky=ky+incy
   90          continue
          go to 140
  100     continue
          sh11=sparam(2)
          sh22=sparam(5)
               do 110 i=1,n
               w=sx(kx)
               z=sy(ky)
               sx(kx)=w*sh11+z
               sy(ky)=-w+sh22*z
               kx=kx+incx
               ky=ky+incy
  110          continue
          go to 140
  120     continue
          sh11=sparam(2)
          sh12=sparam(4)
          sh21=sparam(3)
          sh22=sparam(5)
               do 130 i=1,n
               w=sx(kx)
               z=sy(ky)
               sx(kx)=w*sh11+z*sh12
               sy(ky)=w*sh21+z*sh22
               kx=kx+incx
               ky=ky+incy
  130          continue
  140     continue
          return
          end
