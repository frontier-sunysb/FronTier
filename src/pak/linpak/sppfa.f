      subroutine sppfa(ap,n,info)
      integer n,info
      real ap(1)
c
c     sppfa factors a real symmetric positive definite matrix
c     stored in packed form.
c
c     sppfa is usually called by sppco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for sppco) = (1 + 18/n)*(time for sppfa) .
c
c     on entry
c
c        ap      real (n*(n+1)/2)
c                the packed form of a symmetric matrix  a .  the
c                columns of the upper triangle are stored sequentially
c                in a one-dimensional array of length  n*(n+1)/2 .
c                see comments below for details.
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        ap      an upper triangular matrix  r , stored in packed
c                form, so that  a = trans(r)*r .
c
c        info    integer
c                = 0  for normal return.
c                = k  if the leading minor of order  k  is not
c                     positive definite.
c
c
c     packed storage
c
c          the following program segment will pack the upper
c          triangle of a symmetric matrix.
c
c                k = 0
c                do 20 j = 1, n
c                   do 10 i = 1, j
c                      k = k + 1
c                      ap(k) = a(i,j)
c             10    continue
c             20 continue
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas sdot
c     fortran sqrt
c
c     internal variables
c
      real sdot,t
      real s
      integer j,jj,jm1,k,kj,kk
c     begin block with ...exits to 40
c
c
         jj = 0
         do 30 j = 1, n
            info = j
            s = 0.0e0
            jm1 = j - 1
            kj = jj
            kk = 0
            if (jm1 .lt. 1) go to 20
            do 10 k = 1, jm1
               kj = kj + 1
               t = ap(kj) - sdot(k-1,ap(kk+1),1,ap(jj+1),1)
               kk = kk + k
               t = t/ap(kk)
               ap(kj) = t
               s = s + t*t
   10       continue
   20       continue
            jj = jj + j
            s = ap(jj) - s
c     ......exit
            if (s .le. 0.0e0) go to 40
            ap(jj) = sqrt(s)
   30    continue
         info = 0
   40 continue
      return
      end
