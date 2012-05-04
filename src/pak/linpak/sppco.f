      subroutine sppco(ap,n,rcond,z,info)
      integer n,info
      real ap(1),z(1)
      real rcond
c
c     sppco factors a real symmetric positive definite matrix
c     stored in packed form
c     and estimates the condition of the matrix.
c
c     if  rcond  is not needed, sppfa is slightly faster.
c     to solve  a*x = b , follow sppco by sppsl.
c     to compute  inverse(a)*c , follow sppco by sppsl.
c     to compute  determinant(a) , follow sppco by sppdi.
c     to compute  inverse(a) , follow sppco by sppdi.
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
c                if  info .ne. 0 , the factorization is not complete.
c
c        rcond   real
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.  if info .ne. 0 , rcond is unchanged.
c
c        z       real(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is singular to working precision, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c                if  info .ne. 0 , z  is unchanged.
c
c        info    integer
c                = 0  for normal return.
c                = k  signals an error condition.  the leading minor
c                     of order  k  is not positive definite.
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
c     linpack sppfa
c     blas saxpy,sdot,sscal,sasum
c     fortran abs,amax1,real,sign
c
c     internal variables
c
      real sdot,ek,t,wk,wkm
      real anorm,s,sasum,sm,ynorm
      integer i,ij,j,jm1,j1,k,kb,kj,kk,kp1
c
c
c     find norm of a
c
      j1 = 1
      do 30 j = 1, n
         z(j) = sasum(j,ap(j1),1)
         ij = j1
         j1 = j1 + j
         jm1 = j - 1
         if (jm1 .lt. 1) go to 20
         do 10 i = 1, jm1
            z(i) = z(i) + abs(ap(ij))
            ij = ij + 1
   10    continue
   20    continue
   30 continue
      anorm = 0.0e0
      do 40 j = 1, n
         anorm = amax1(anorm,z(j))
   40 continue
c
c     factor
c
      call sppfa(ap,n,info)
      if (info .ne. 0) go to 180
c
c        rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c        estimate = norm(z)/norm(y) where  a*z = y  and  a*y = e .
c        the components of  e  are chosen to cause maximum local
c        growth in the elements of w  where  trans(r)*w = e .
c        the vectors are frequently rescaled to avoid overflow.
c
c        solve trans(r)*w = e
c
         ek = 1.0e0
         do 50 j = 1, n
            z(j) = 0.0e0
   50    continue
         kk = 0
         do 110 k = 1, n
            kk = kk + k
            if (z(k) .ne. 0.0e0) ek = sign(ek,-z(k))
            if (abs(ek-z(k)) .le. ap(kk)) go to 60
               s = ap(kk)/abs(ek-z(k))
               call sscal(n,s,z,1)
               ek = s*ek
   60       continue
            wk = ek - z(k)
            wkm = -ek - z(k)
            s = abs(wk)
            sm = abs(wkm)
            wk = wk/ap(kk)
            wkm = wkm/ap(kk)
            kp1 = k + 1
            kj = kk + k
            if (kp1 .gt. n) go to 100
               do 70 j = kp1, n
                  sm = sm + abs(z(j)+wkm*ap(kj))
                  z(j) = z(j) + wk*ap(kj)
                  s = s + abs(z(j))
                  kj = kj + j
   70          continue
               if (s .ge. sm) go to 90
                  t = wkm - wk
                  wk = wkm
                  kj = kk + k
                  do 80 j = kp1, n
                     z(j) = z(j) + t*ap(kj)
                     kj = kj + j
   80             continue
   90          continue
  100       continue
            z(k) = wk
  110    continue
         s = 1.0e0/sasum(n,z,1)
         call sscal(n,s,z,1)
c
c        solve r*y = w
c
         do 130 kb = 1, n
            k = n + 1 - kb
            if (abs(z(k)) .le. ap(kk)) go to 120
               s = ap(kk)/abs(z(k))
               call sscal(n,s,z,1)
  120       continue
            z(k) = z(k)/ap(kk)
            kk = kk - k
            t = -z(k)
            call saxpy(k-1,t,ap(kk+1),1,z(1),1)
  130    continue
         s = 1.0e0/sasum(n,z,1)
         call sscal(n,s,z,1)
c
         ynorm = 1.0e0
c
c        solve trans(r)*v = y
c
         do 150 k = 1, n
            z(k) = z(k) - sdot(k-1,ap(kk+1),1,z(1),1)
            kk = kk + k
            if (abs(z(k)) .le. ap(kk)) go to 140
               s = ap(kk)/abs(z(k))
               call sscal(n,s,z,1)
               ynorm = s*ynorm
  140       continue
            z(k) = z(k)/ap(kk)
  150    continue
         s = 1.0e0/sasum(n,z,1)
         call sscal(n,s,z,1)
         ynorm = s*ynorm
c
c        solve r*z = v
c
         do 170 kb = 1, n
            k = n + 1 - kb
            if (abs(z(k)) .le. ap(kk)) go to 160
               s = ap(kk)/abs(z(k))
               call sscal(n,s,z,1)
               ynorm = s*ynorm
  160       continue
            z(k) = z(k)/ap(kk)
            kk = kk - k
            t = -z(k)
            call saxpy(k-1,t,ap(kk+1),1,z(1),1)
  170    continue
c        make znorm = 1.0
         s = 1.0e0/sasum(n,z,1)
         call sscal(n,s,z,1)
         ynorm = s*ynorm
c
         if (anorm .ne. 0.0e0) rcond = ynorm/anorm
         if (anorm .eq. 0.0e0) rcond = 0.0e0
  180 continue
      return
      end
