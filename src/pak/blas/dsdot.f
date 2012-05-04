      double precision function dsdot(n,sx,incx,sy,incy)
c
c     returns d.p. dot product accumulated in d.p., for s.p. sx and sy
c     dsdot = sum for i = 0 to n-1 of  sx(lx+i*incx) * sy(ly+i*incy),
c     where lx = 1 if incx .ge. 0, else lx = (-incx)*n, and ly is
c     defined in a similar way using incy.
c
      real sx(1),sy(1)
c
      dsdot = 0.d0
      if(n .le. 0)return
      if(incx.eq.incy.and.incx.gt.0) go to 20
      kx = 1
      ky = 1
      if(incx.lt.0) kx = 1+(1-n)*incx
      if(incy.lt.0) ky = 1+(1-n)*incy
          do 10 i = 1,n
          dsdot = dsdot + dble(sx(kx))*dble(sy(ky))
          kx = kx + incx
          ky = ky + incy
   10 continue
      return
   20 continue
      ns = n*incx
          do 30 i=1,ns,incx
          dsdot = dsdot + dble(sx(i))*dble(sy(i))
   30     continue
      return
      end
