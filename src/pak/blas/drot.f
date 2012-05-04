      subroutine drot(n,dx,incx,dy,incy,dc,ds)
c
c     multiply the 2 x 2 matrix  ( dc ds) times the 2 x n matrix (dx**t)
c                                (-ds dc)                        (dy**t)
c     where **t indicates transpose.    the elements of dx are in
c     dx(lx+i*incx), i = 0 to n-1, where lx = 1 if incx .ge. 0, else
c     lx = (-incx)*n, and similarly for dy using ly and incy.
      double precision dx,dy,dc,ds,zero,one,w,z
      dimension dx(1),dy(1)
c
      data zero,one/0.d0,1.d0/
      if(n .le. 0 .or. (ds .eq. zero .and. dc .eq. one)) go to 40
      if(.not. (incx .eq. incy .and. incx .gt. 0)) go to 20
c
           nsteps=incx*n
           do 10 i=1,nsteps,incx
                w=dx(i)
                z=dy(i)
                dx(i)=dc*w+ds*z
                dy(i)=-ds*w+dc*z
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
                w=dx(kx)
                z=dy(ky)
                dx(kx)=dc*w+ds*z
                dy(ky)=-ds*w+dc*z
                kx=kx+incx
                ky=ky+incy
   30           continue
   40 continue
c
      return
      end
