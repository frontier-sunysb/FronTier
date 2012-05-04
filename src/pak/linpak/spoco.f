      subroutine spoco(a,lda,n,rcond,z,info)
      integer lda,n,info
      real a(lda,1),z(1)
      real rcond
c
c     spoco factors a real symmetric positive definite matrix
c     and estimates the condition of the matrix.
c
c     if  rcond  is not needed, spofa is slightly faster.
c     to solve  a*x = b , follow spoco by sposl.
c     to compute  inverse(a)*c , follow spoco by sposl.
c     to compute  determinant(a) , follow spoco by spodi.
c     to compute  inverse(a) , follow spoco by spodi.
c
c     on entry
c
c        a       real(lda, n)
c                the symmetric matrix to be factored.  only the
c                diagonal and upper triangle are used.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix  r  so that  a = trans(r)*r
c                where  trans(r)  is the transpose.
c                the strict lower triangle is unaltered.
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
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c                if  info .ne. 0 , z  is unchanged.
c
c        info    integer
c                = 0  for normal return.
c                = k  signals an error condition.  the leading minor
c                     of order  k  is not positive definite.
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack spofa
c     blas saxpy,sdot,sscal,sasum
c     fortran abs,amax1,real,sign
c
c     internal variables
c
      real sdot,ek,t,wk,wkm
      real anorm,s,sasum,sm,ynorm
      integer i,j,jm1,k,kb,kp1
c
c
c     find norm of a using only upper half
c
      do 30 j = 1, n
         z(j) = sasum(j,a(1,j),1)
         jm1 = j - 1
         if (jm1 .lt. 1) go to 20
         do 10 i = 1, jm1
            z(i) = z(i) + abs(a(i,j))
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
      call spofa(a,lda,n,info)
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
         do 110 k = 1, n
            if (z(k) .ne. 0.0e0) ek = sign(ek,-z(k))
            if (abs(ek-z(k)) .le. a(k,k)) go to 60
               s = a(k,k)/abs(ek-z(k))
               call sscal(n,s,z,1)
               ek = s*ek
   60       continue
            wk = ek - z(k)
            wkm = -ek - z(k)
            s = abs(wk)
            sm = abs(wkm)
            wk = wk/a(k,k)
            wkm = wkm/a(k,k)
            kp1 = k + 1
            if (kp1 .gt. n) go to 100
               do 70 j = kp1, n
                  sm = sm + abs(z(j)+wkm*a(k,j))
                  z(j) = z(j) + wk*a(k,j)
                  s = s + abs(z(j))
   70          continue
               if (s .ge. sm) go to 90
                  t = wkm - wk
                  wk = wkm
                  do 80 j = kp1, n
                     z(j) = z(j) + t*a(k,j)
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
            if (abs(z(k)) .le. a(k,k)) go to 120
               s = a(k,k)/abs(z(k))
               call sscal(n,s,z,1)
  120       continue
            z(k) = z(k)/a(k,k)
            t = -z(k)
            call saxpy(k-1,t,a(1,k),1,z(1),1)
  130    continue
         s = 1.0e0/sasum(n,z,1)
         call sscal(n,s,z,1)
c
         ynorm = 1.0e0
c
c        solve trans(r)*v = y
c
         do 150 k = 1, n
            z(k) = z(k) - sdot(k-1,a(1,k),1,z(1),1)
            if (abs(z(k)) .le. a(k,k)) go to 140
               s = a(k,k)/abs(z(k))
               call sscal(n,s,z,1)
               ynorm = s*ynorm
  140       continue
            z(k) = z(k)/a(k,k)
  150    continue
         s = 1.0e0/sasum(n,z,1)
         call sscal(n,s,z,1)
         ynorm = s*ynorm
c
c        solve r*z = v
c
         do 170 kb = 1, n
            k = n + 1 - kb
            if (abs(z(k)) .le. a(k,k)) go to 160
               s = a(k,k)/abs(z(k))
               call sscal(n,s,z,1)
               ynorm = s*ynorm
  160       continue
            z(k) = z(k)/a(k,k)
            t = -z(k)
            call saxpy(k-1,t,a(1,k),1,z(1),1)
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
