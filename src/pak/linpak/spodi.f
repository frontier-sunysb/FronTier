      subroutine spodi(a,lda,n,det,job)
      integer lda,n,job
      real a(lda,1)
      real det(2)
c
c     spodi computes the determinant and inverse of a certain
c     real symmetric positive definite matrix (see below)
c     using the factors computed by spoco, spofa or sqrdc.
c
c     on entry
c
c        a       real(lda, n)
c                the output  a  from spoco or spofa
c                or the output  x  from sqrdc.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        a       if spoco or spofa was used to factor  a  then
c                spodi produces the upper half of inverse(a) .
c                if sqrdc was used to decompose  x  then
c                spodi produces the upper half of inverse(trans(x)*x)
c                where trans(x) is the transpose.
c                elements of  a  below the diagonal are unchanged.
c                if the units digit of job is zero,  a  is unchanged.
c
c        det     real(2)
c                determinant of  a  or of  trans(x)*x  if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. det(1) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if spoco or spofa has set info .eq. 0 .
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas saxpy,sscal
c     fortran mod
c
c     internal variables
c
      real t
      real s
      integer i,j,jm1,k,kp1
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
         det(1) = 1.0e0
         det(2) = 0.0e0
         s = 10.0e0
         do 50 i = 1, n
            det(1) = a(i,i)**2*det(1)
c        ...exit
            if (det(1) .eq. 0.0e0) go to 60
   10       if (det(1) .ge. 1.0e0) go to 20
               det(1) = s*det(1)
               det(2) = det(2) - 1.0e0
            go to 10
   20       continue
   30       if (det(1) .lt. s) go to 40
               det(1) = det(1)/s
               det(2) = det(2) + 1.0e0
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(r)
c
      if (mod(job,10) .eq. 0) go to 140
         do 100 k = 1, n
            a(k,k) = 1.0e0/a(k,k)
            t = -a(k,k)
            call sscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = a(k,j)
               a(k,j) = 0.0e0
               call saxpy(k,t,a(1,k),1,a(1,j),1)
   80       continue
   90       continue
  100    continue
c
c        form  inverse(r) * trans(inverse(r))
c
         do 130 j = 1, n
            jm1 = j - 1
            if (jm1 .lt. 1) go to 120
            do 110 k = 1, jm1
               t = a(k,j)
               call saxpy(k,t,a(1,j),1,a(1,k),1)
  110       continue
  120       continue
            t = a(j,j)
            call sscal(j,t,a(1,j),1)
  130    continue
  140 continue
      return
      end
