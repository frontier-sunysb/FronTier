      subroutine ctrsl(t,ldt,n,b,job,info)
      integer ldt,n,job,info
      complex t(ldt,1),b(1)
c
c
c     ctrsl solves systems of the form
c
c                   t * x = b
c     or
c                   ctrans(t) * x = b
c
c     where t is a triangular matrix of order n. here ctrans(t)
c     denotes the conjugate transpose of the matrix t.
c
c     on entry
c
c         t         complex(ldt,n)
c                   t contains the matrix of the system. the zero
c                   elements of the matrix are not referenced, and
c                   the corresponding elements of the array can be
c                   used to store other information.
c
c         ldt       integer
c                   ldt is the leading dimension of the array t.
c
c         n         integer
c                   n is the order of the system.
c
c         b         complex(n).
c                   b contains the right hand side of the system.
c
c         job       integer
c                   job specifies what kind of system is to be solved.
c                   if job is
c
c                        00   solve t*x=b, t lower triangular,
c                        01   solve t*x=b, t upper triangular,
c                        10   solve ctrans(t)*x=b, t lower triangular,
c                        11   solve ctrans(t)*x=b, t upper triangular.
c
c     on return
c
c         b         b contains the solution, if info .eq. 0.
c                   otherwise b is unaltered.
c
c         info      integer
c                   info contains zero if the system is nonsingular.
c                   otherwise info contains the index of
c                   the first zero diagonal element of t.
c
c     linpack. this version dated 08/14/78 .
c     g. w. stewart, university of maryland, argonne national lab.
c
c     subroutines and functions
c
c     blas caxpy,cdotc
c     fortran abs,aimag,conjg,mod,real
c
c     internal variables
c
      complex cdotc,temp
      integer case,j,jj
c
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c     begin block permitting ...exits to 150
c
c        check for zero diagonal elements.
c
         do 10 info = 1, n
c     ......exit
            if (cabs1(t(info,info)) .eq. 0.0e0) go to 150
   10    continue
         info = 0
c
c        determine the task and go to it.
c
         case = 1
         if (mod(job,10) .ne. 0) case = 2
         if (mod(job,100)/10 .ne. 0) case = case + 2
         go to (20,50,80,110), case
c
c        solve t*x=b for t lower triangular
c
   20    continue
            b(1) = b(1)/t(1,1)
            if (n .lt. 2) go to 40
            do 30 j = 2, n
               temp = -b(j-1)
               call caxpy(n-j+1,temp,t(j,j-1),1,b(j),1)
               b(j) = b(j)/t(j,j)
   30       continue
   40       continue
         go to 140
c
c        solve t*x=b for t upper triangular.
c
   50    continue
            b(n) = b(n)/t(n,n)
            if (n .lt. 2) go to 70
            do 60 jj = 2, n
               j = n - jj + 1
               temp = -b(j+1)
               call caxpy(j,temp,t(1,j+1),1,b(1),1)
               b(j) = b(j)/t(j,j)
   60       continue
   70       continue
         go to 140
c
c        solve ctrans(t)*x=b for t lower triangular.
c
   80    continue
            b(n) = b(n)/conjg(t(n,n))
            if (n .lt. 2) go to 100
            do 90 jj = 2, n
               j = n - jj + 1
               b(j) = b(j) - cdotc(jj-1,t(j+1,j),1,b(j+1),1)
               b(j) = b(j)/conjg(t(j,j))
   90       continue
  100       continue
         go to 140
c
c        solve ctrans(t)*x=b for t upper triangular.
c
  110    continue
            b(1) = b(1)/conjg(t(1,1))
            if (n .lt. 2) go to 130
            do 120 j = 2, n
               b(j) = b(j) - cdotc(j-1,t(1,j),1,b(1),1)
               b(j) = b(j)/conjg(t(j,j))
  120       continue
  130       continue
  140    continue
  150 continue
      return
      end
