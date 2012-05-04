      subroutine ctrco(t,ldt,n,rcond,z,job)
      integer ldt,n,job
      complex t(ldt,1),z(1)
      real rcond
c
c     ctrco estimates the condition of a complex triangular matrix.
c
c     on entry
c
c        t       complex(ldt,n)
c                t contains the triangular matrix. the zero
c                elements of the matrix are not referenced, and
c                the corresponding elements of the array can be
c                used to store other information.
c
c        ldt     integer
c                ldt is the leading dimension of the array t.
c
c        n       integer
c                n is the order of the system.
c
c        job     integer
c                = 0         t  is lower triangular.
c                = nonzero   t  is upper triangular.
c
c     on return
c
c        rcond   real
c                an estimate of the reciprocal condition of  t .
c                for the system  t*x = b , relative perturbations
c                in  t  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  t  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       complex(n)
c                a work vector whose contents are usually unimportant.
c                if  t  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas caxpy,csscal,scasum
c     fortran abs,aimag,amax1,cmplx,conjg,real
c
c     internal variables
c
      complex w,wk,wkm,ek
      real tnorm,ynorm,s,sm,scasum
      integer i1,j,j1,j2,k,kk,l
      logical lower
      complex zdum,zdum1,zdum2,csign1
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
      csign1(zdum1,zdum2) = cabs1(zdum1)*(zdum2/cabs1(zdum2))
c
      lower = job .eq. 0
c
c     compute 1-norm of t
c
      tnorm = 0.0e0
      do 10 j = 1, n
         l = j
         if (lower) l = n + 1 - j
         i1 = 1
         if (lower) i1 = j
         tnorm = amax1(tnorm,scasum(l,t(i1,j),1))
   10 continue
c
c     rcond = 1/(norm(t)*(estimate of norm(inverse(t)))) .
c     estimate = norm(z)/norm(y) where  t*z = y  and  ctrans(t)*y = e .
c     ctrans(t)  is the conjugate transpose of t .
c     the components of  e  are chosen to cause maximum local
c     growth in the elements of y .
c     the vectors are frequently rescaled to avoid overflow.
c
c     solve ctrans(t)*y = e
c
      ek = (1.0e0,0.0e0)
      do 20 j = 1, n
         z(j) = (0.0e0,0.0e0)
   20 continue
      do 100 kk = 1, n
         k = kk
         if (lower) k = n + 1 - kk
         if (cabs1(z(k)) .ne. 0.0e0) ek = csign1(ek,-z(k))
         if (cabs1(ek-z(k)) .le. cabs1(t(k,k))) go to 30
            s = cabs1(t(k,k))/cabs1(ek-z(k))
            call csscal(n,s,z,1)
            ek = cmplx(s,0.0e0)*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = cabs1(wk)
         sm = cabs1(wkm)
         if (cabs1(t(k,k)) .eq. 0.0e0) go to 40
            wk = wk/conjg(t(k,k))
            wkm = wkm/conjg(t(k,k))
         go to 50
   40    continue
            wk = (1.0e0,0.0e0)
            wkm = (1.0e0,0.0e0)
   50    continue
         if (kk .eq. n) go to 90
            j1 = k + 1
            if (lower) j1 = 1
            j2 = n
            if (lower) j2 = k - 1
            do 60 j = j1, j2
               sm = sm + cabs1(z(j)+wkm*conjg(t(k,j)))
               z(j) = z(j) + wk*conjg(t(k,j))
               s = s + cabs1(z(j))
   60       continue
            if (s .ge. sm) go to 80
               w = wkm - wk
               wk = wkm
               do 70 j = j1, j2
                  z(j) = z(j) + w*conjg(t(k,j))
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
c
      ynorm = 1.0e0
c
c     solve t*z = y
c
      do 130 kk = 1, n
         k = n + 1 - kk
         if (lower) k = kk
         if (cabs1(z(k)) .le. cabs1(t(k,k))) go to 110
            s = cabs1(t(k,k))/cabs1(z(k))
            call csscal(n,s,z,1)
            ynorm = s*ynorm
  110    continue
         if (cabs1(t(k,k)) .ne. 0.0e0) z(k) = z(k)/t(k,k)
         if (cabs1(t(k,k)) .eq. 0.0e0) z(k) = (1.0e0,0.0e0)
         i1 = 1
         if (lower) i1 = k + 1
         if (kk .ge. n) go to 120
            w = -z(k)
            call caxpy(n-kk,w,t(i1,k),1,z(i1),1)
  120    continue
  130 continue
c     make znorm = 1.0
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (tnorm .ne. 0.0e0) rcond = ynorm/tnorm
      if (tnorm .eq. 0.0e0) rcond = 0.0e0
      return
      end
