      subroutine cpoco(a,lda,n,rcond,z,info)
      integer lda,n,info
      complex a(lda,1),z(1)
      real rcond
c
c     cpoco factors a complex hermitian positive definite matrix
c     and estimates the condition of the matrix.
c
c     if  rcond  is not needed, cpofa is slightly faster.
c     to solve  a*x = b , follow cpoco by cposl.
c     to compute  inverse(a)*c , follow cpoco by cposl.
c     to compute  determinant(a) , follow cpoco by cpodi.
c     to compute  inverse(a) , follow cpoco by cpodi.
c
c     on entry
c
c        a       complex(lda, n)
c                the hermitian matrix to be factored.  only the
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
c        a       an upper triangular matrix  r  so that  a =
c                ctrans(r)*r where  ctrans(r)  is the conjugate
c                transpose.  the strict lower triangle is unaltered.
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
c        z       complex(n)
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
c     linpack cpofa
c     blas caxpy,cdotc,csscal,scasum
c     fortran abs,aimag,amax1,cmplx,conjg,real
c
c     internal variables
c
      complex cdotc,ek,t,wk,wkm
      real anorm,s,scasum,sm,ynorm
      integer i,j,jm1,k,kb,kp1
c
      complex zdum,zdum2,csign1
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
      csign1(zdum,zdum2) = cabs1(zdum)*(zdum2/cabs1(zdum2))
c
c     find norm of a using only upper half
c
      do 30 j = 1, n
         z(j) = cmplx(scasum(j,a(1,j),1),0.0e0)
         jm1 = j - 1
         if (jm1 .lt. 1) go to 20
         do 10 i = 1, jm1
            z(i) = cmplx(real(z(i))+cabs1(a(i,j)),0.0e0)
   10    continue
   20    continue
   30 continue
      anorm = 0.0e0
      do 40 j = 1, n
         anorm = amax1(anorm,real(z(j)))
   40 continue
c
c     factor
c
      call cpofa(a,lda,n,info)
      if (info .ne. 0) go to 180
c
c        rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c        estimate = norm(z)/norm(y) where  a*z = y  and  a*y = e .
c        the components of  e  are chosen to cause maximum local
c        growth in the elements of w  where  ctrans(r)*w = e .
c        the vectors are frequently rescaled to avoid overflow.
c
c        solve ctrans(r)*w = e
c
         ek = (1.0e0,0.0e0)
         do 50 j = 1, n
            z(j) = (0.0e0,0.0e0)
   50    continue
         do 110 k = 1, n
            if (cabs1(z(k)) .ne. 0.0e0) ek = csign1(ek,-z(k))
            if (cabs1(ek-z(k)) .le. real(a(k,k))) go to 60
               s = real(a(k,k))/cabs1(ek-z(k))
               call csscal(n,s,z,1)
               ek = cmplx(s,0.0e0)*ek
   60       continue
            wk = ek - z(k)
            wkm = -ek - z(k)
            s = cabs1(wk)
            sm = cabs1(wkm)
            wk = wk/a(k,k)
            wkm = wkm/a(k,k)
            kp1 = k + 1
            if (kp1 .gt. n) go to 100
               do 70 j = kp1, n
                  sm = sm + cabs1(z(j)+wkm*conjg(a(k,j)))
                  z(j) = z(j) + wk*conjg(a(k,j))
                  s = s + cabs1(z(j))
   70          continue
               if (s .ge. sm) go to 90
                  t = wkm - wk
                  wk = wkm
                  do 80 j = kp1, n
                     z(j) = z(j) + t*conjg(a(k,j))
   80             continue
   90          continue
  100       continue
            z(k) = wk
  110    continue
         s = 1.0e0/scasum(n,z,1)
         call csscal(n,s,z,1)
c
c        solve r*y = w
c
         do 130 kb = 1, n
            k = n + 1 - kb
            if (cabs1(z(k)) .le. real(a(k,k))) go to 120
               s = real(a(k,k))/cabs1(z(k))
               call csscal(n,s,z,1)
  120       continue
            z(k) = z(k)/a(k,k)
            t = -z(k)
            call caxpy(k-1,t,a(1,k),1,z(1),1)
  130    continue
         s = 1.0e0/scasum(n,z,1)
         call csscal(n,s,z,1)
c
         ynorm = 1.0e0
c
c        solve ctrans(r)*v = y
c
         do 150 k = 1, n
            z(k) = z(k) - cdotc(k-1,a(1,k),1,z(1),1)
            if (cabs1(z(k)) .le. real(a(k,k))) go to 140
               s = real(a(k,k))/cabs1(z(k))
               call csscal(n,s,z,1)
               ynorm = s*ynorm
  140       continue
            z(k) = z(k)/a(k,k)
  150    continue
         s = 1.0e0/scasum(n,z,1)
         call csscal(n,s,z,1)
         ynorm = s*ynorm
c
c        solve r*z = v
c
         do 170 kb = 1, n
            k = n + 1 - kb
            if (cabs1(z(k)) .le. real(a(k,k))) go to 160
               s = real(a(k,k))/cabs1(z(k))
               call csscal(n,s,z,1)
               ynorm = s*ynorm
  160       continue
            z(k) = z(k)/a(k,k)
            t = -z(k)
            call caxpy(k-1,t,a(1,k),1,z(1),1)
  170    continue
c        make znorm = 1.0
         s = 1.0e0/scasum(n,z,1)
         call csscal(n,s,z,1)
         ynorm = s*ynorm
c
         if (anorm .ne. 0.0e0) rcond = ynorm/anorm
         if (anorm .eq. 0.0e0) rcond = 0.0e0
  180 continue
      return
      end
