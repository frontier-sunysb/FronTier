      subroutine cchud(r,ldr,p,x,z,ldz,nz,y,rho,c,s)
      integer ldr,p,ldz,nz
      real rho(1),c(1)
      complex r(ldr,1),x(1),z(ldz,1),y(1),s(1)
c
c     cchud updates an augmented cholesky decomposition of the
c     triangular part of an augmented qr decomposition.  specifically,
c     given an upper triangular matrix r of order p, a row vector
c     x, a column vector z, and a scalar y, cchud determines a
c     untiary matrix u and a scalar zeta such that
c
c
c                              (r  z)     (rr   zz )
c                         u  * (    )  =  (        ) ,
c                              (x  y)     ( 0  zeta)
c
c     where rr is upper triangular.  if r and z have been
c     obtained from the factorization of a least squares
c     problem, then rr and zz are the factors corresponding to
c     the problem with the observation (x,y) appended.  in this
c     case, if rho is the norm of the residual vector, then the
c     norm of the residual vector of the updated problem is
c     sqrt(rho**2 + zeta**2).  cchud will simultaneously update
c     several triplets (z,y,rho).
c     for a less terse description of what cchud does and how
c     it may be applied see the linpack guide.
c
c     the matrix u is determined as the product u(p)*...*u(1),
c     where u(i) is a rotation in the (i,p+1) plane of the
c     form
c
c                       (     c(i)      s(i) )
c                       (                    ) .
c                       ( -conjg(s(i))  c(i) )
c
c     the rotations are chosen so that c(i) is real.
c
c     on entry
c
c         r      complex(ldr,p), where ldr .ge. p.
c                r contains the upper triangular matrix
c                that is to be updated.  the part of r
c                below the diagonal is not referenced.
c
c         ldr    integer.
c                ldr is the leading dimension of the array r.
c
c         p      integer.
c                p is the order of the matrix r.
c
c         x      complex(p).
c                x contains the row to be added to r.  x is
c                not altered by cchud.
c
c         z      complex(ldz,nz), where ldz .ge. p.
c                z is an array containing nz p-vectors to
c                be updated with r.
c
c         ldz    integer.
c                ldz is the leading dimension of the array z.
c
c         nz     integer.
c                nz is the number of vectors to be updated
c                nz may be zero, in which case z, y, and rho
c                are not referenced.
c
c         y      complex(nz).
c                y contains the scalars for updating the vectors
c                z.  y is not altered by cchud.
c
c         rho    real(nz).
c                rho contains the norms of the residual
c                vectors that are to be updated.  if rho(j)
c                is negative, it is left unaltered.
c
c     on return
c
c         rc
c         rho    contain the updated quantities.
c         z
c
c         c      real(p).
c                c contains the cosines of the transforming
c                rotations.
c
c         s      complex(p).
c                s contains the sines of the transforming
c                rotations.
c
c     linpack. this version dated 08/14/78 .
c     g.w. stewart, university of maryland, argonne national lab.
c
c     cchud uses the following functions and subroutines.
c
c     extended blas crotg
c     fortran conjg,sqrt
c
      integer i,j,jm1
      real azeta,scale
      complex t,xj,zeta
c
c     update r.
c
      do 30 j = 1, p
         xj = x(j)
c
c        apply the previous rotations.
c
         jm1 = j - 1
         if (jm1 .lt. 1) go to 20
         do 10 i = 1, jm1
            t = c(i)*r(i,j) + s(i)*xj
            xj = c(i)*xj - conjg(s(i))*r(i,j)
            r(i,j) = t
   10    continue
   20    continue
c
c        compute the next rotation.
c
         call crotg(r(j,j),xj,c(j),s(j))
   30 continue
c
c     if required, update z and rho.
c
      if (nz .lt. 1) go to 70
      do 60 j = 1, nz
         zeta = y(j)
         do 40 i = 1, p
            t = c(i)*z(i,j) + s(i)*zeta
            zeta = c(i)*zeta - conjg(s(i))*z(i,j)
            z(i,j) = t
   40    continue
         azeta = cabs(zeta)
         if (azeta .eq. 0.0e0 .or. rho(j) .lt. 0.0e0) go to 50
            scale = azeta + rho(j)
            rho(j) = scale*sqrt((azeta/scale)**2+(rho(j)/scale)**2)
   50    continue
   60 continue
   70 continue
      return
      end
