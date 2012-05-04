      real function sdsdot(n,sb,sx,incx,sy,incy)
c
c     returns s.p. result with dot product accumulated in d.p.
c     sdsdot = sb + sum for i = 0 to n-1 of sx(lx+i*incx)*sy(ly+i*incy),
c     where lx = 1 if incx .ge. 0, else lx = (-incx)*n, and ly is
c     defined in a similar way using incy.
c
      real              sx(1),sy(1),sb
      double precision dsdot
c
      dsdot = dble(sb)
      if(n .le. 0) go to 30
      if(incx.eq.incy.and.incx.gt.0) go to 40
      kx = 1
      ky = 1
      if(incx.lt.0) kx = 1+(1-n)*incx
      if(incy.lt.0) ky = 1+(1-n)*incy
          do 10 i = 1,n
          dsdot = dsdot + dble(sx(kx))*dble(sy(ky))
          kx = kx + incx
          ky = ky + incy
   10     continue
   30 sdsdot = sngl(dsdot)
      return
   40 continue
      ns = n*incx
          do 50 i=1,ns,incx
          dsdot = dsdot + dble(sx(i))*dble(sy(i))
   50     continue
      sdsdot = sngl(dsdot)
      return
      end
