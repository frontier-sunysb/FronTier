      subroutine srot(n,sx,incx,sy,incy,sc,ss)
c
c     multiply the 2 x 2 matrix  ( sc ss) times the 2 x n matrix (sx**t)
c                                (-ss sc)                        (sy**t)
c     where **t indicates transpose.    the elements of sx are in
c     sx(lx+i*incx), i = 0 to n-1, where lx = 1 if incx .ge. 0, else
c     lx = (-incx)*n, and similarly for sy using ly and incy.
      real             sx,sy,sc,ss,zero,one,w,z
      dimension sx(1),sy(1)
c
      data zero,one/0.e0,1.e0/
      if(n .le. 0 .or. (ss .eq. zero .and. sc .eq. one)) go to 40
      if(.not. (incx .eq. incy .and. incx .gt. 0)) go to 20
c
           nsteps=incx*n
           do 10 i=1,nsteps,incx
                w=sx(i)
                z=sy(i)
                sx(i)=sc*w+ss*z
                sy(i)=-ss*w+sc*z
   10           continue
           go to 40
c
   20 continue
           kx=1
           ky=1
c
           if(incx .lt. 0) kx=1-(n-1)*incx
           if(incy .lt. 0) ky=1-(n-1)*incy
c
           do 30 i=1,n
                w=sx(kx)
                z=sy(ky)
                sx(kx)=sc*w+ss*z
                sy(ky)=-ss*w+sc*z
                kx=kx+incx
                ky=ky+incy
   30           continue
   40 continue
c
      return
      end
