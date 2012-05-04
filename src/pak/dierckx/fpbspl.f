      subroutine fpbspl(t,n,k,x,l,h)
      implicit none
c  subroutine fpbspl evaluates the (k+1) non-zero b-splines of
c  degree k at t(l) <= x < t(l+1) using the stable recurrence
c  relation of de boor and cox.
c  ..
c  ..scalar arguments..
      double precision x
      integer n,k,l
c  ..array arguments..
      double precision t(n),h(6)
c  ..local scalars..
      double precision f,one
      integer i,j,li,lj
c  ..local arrays..
      double precision hh(5)
c  ..
      one = 0.1d+01
      h(1) = one
      do j=1,k
        do i=1,j
          hh(i) = h(i)
        enddo
        h(1) = 0.
        do i=1,j
          li = l+i
          lj = li-j
          f = hh(i)/(t(li)-t(lj))
          h(i) = h(i)+f*(t(li)-x)
          h(i+1) = f*(x-t(lj))
        enddo
      enddo
      return
      end
