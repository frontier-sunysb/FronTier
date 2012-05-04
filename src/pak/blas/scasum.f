      function scasum(n,cx,incx)
c     returns sums of magnitudes of real and imaginary parts of
c     components of cx.  note that this is not the l1 norm of cx.
c     casum = sum from 0 to n-1 of abs(real(cx(1+i*incx))) +
c             abs(imag(cx(1+i*incx)))
c
      complex cx(1)
c
      scasum=0.
      if(n .le. 0) return
      ns = n*incx
      do 10 i=1,ns,incx
          scasum = scasum + abs(real(cx(i))) + abs(aimag(cx(i)))
   10     continue
      return
      end
