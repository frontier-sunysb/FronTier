      subroutine srotmg (sd1,sd2,sx1,sy1,sparam)
c
c     construct the modified givens transformation matrix h which zeros
c     the second component of the 2-vector  (sqrt(sd1)*sx1,sqrt(sd2)*
c     sy2)**t.
c     with sparam(1)=sflag, h has one of the following forms..
c
c     sflag=-1.e0     sflag=0.e0        sflag=1.e0     sflag=-2.e0
c
c       (sh11  sh12)    (1.e0  sh12)    (sh11  1.e0)    (1.e0  0.e0)
c     h=(          )    (          )    (          )    (          )
c       (sh21  sh22),   (sh21  1.e0),   (-1.e0 sh22),   (0.e0  1.e0).
c     locations 2-4 of sparam contain sh11,sh21,sh12, and sh22
c     respectively. (values of 1.e0, -1.e0, or 0.e0 implied by the
c     value of sparam(1) are not stored in sparam.)
c
c     the values of gamsq and rgamsq set in the data statement may be
c     inexact.  this is ok as they are only used for testing the size
c     of sd1 and sd2.  all actual scaling of data is done using gam.
c
      dimension sparam(5)
c
      data zero,one,two /0.e0,1.e0,2.e0/
      data gam,gamsq,rgamsq/4096.e0,1.67772e7,5.96046e-8/
      if(.not. sd1 .lt. zero) go to 10
c       go zero-h-d-and-sx1..
          go to 60
   10 continue
c     case-sd1-nonnegative
      sp2=sd2*sy1
      if(.not. sp2 .eq. zero) go to 20
          sflag=-two
          go to 260
c     regular-case..
   20 continue
      sp1=sd1*sx1
      sq2=sp2*sy1
      sq1=sp1*sx1
c
      if(.not. abs(sq1) .gt. abs(sq2)) go to 40
          sh21=-sy1/sx1
          sh12=sp2/sp1
c
          su=one-sh12*sh21
c
          if(.not. su .le. zero) go to 30
c         go zero-h-d-and-sx1..
               go to 60
   30     continue
               sflag=zero
               sd1=sd1/su
               sd2=sd2/su
               sx1=sx1*su
c         go scale-check..
               go to 100
   40 continue
          if(.not. sq2 .lt. zero) go to 50
c         go zero-h-d-and-sx1..
               go to 60
   50     continue
               sflag=one
               sh11=sp1/sp2
               sh22=sx1/sy1
               su=one+sh11*sh22
               stemp=sd2/su
               sd2=sd1/su
               sd1=stemp
               sx1=sy1*su
c         go scale-check
               go to 100
c     procedure..zero-h-d-and-sx1..
   60 continue
          sflag=-one
          sh11=zero
          sh12=zero
          sh21=zero
          sh22=zero
c
          sd1=zero
          sd2=zero
          sx1=zero
c         return..
          go to 220
c     procedure..fix-h..
   70 continue
      if(.not. sflag .ge. zero) go to 90
c
          if(.not. sflag .eq. zero) go to 80
          sh11=one
          sh22=one
          sflag=-one
          go to 90
   80     continue
          sh21=-one
          sh12=one
          sflag=-one
   90 continue
      go to igo,(120,150,180,210)
c     procedure..scale-check
  100 continue
  110     continue
          if(.not. sd1 .le. rgamsq) go to 130
               if(sd1 .eq. zero) go to 160
               assign 120 to igo
c              fix-h..
               go to 70
  120          continue
               sd1=sd1*gam**2
               sx1=sx1/gam
               sh11=sh11/gam
               sh12=sh12/gam
          go to 110
  130 continue
  140     continue
          if(.not. sd1 .ge. gamsq) go to 160
               assign 150 to igo
c              fix-h..
               go to 70
  150          continue
               sd1=sd1/gam**2
               sx1=sx1*gam
               sh11=sh11*gam
               sh12=sh12*gam
          go to 140
  160 continue
  170     continue
          if(.not. abs(sd2) .le. rgamsq) go to 190
               if(sd2 .eq. zero) go to 220
               assign 180 to igo
c              fix-h..
               go to 70
  180          continue
               sd2=sd2*gam**2
               sh21=sh21/gam
               sh22=sh22/gam
          go to 170
  190 continue
  200     continue
          if(.not. abs(sd2) .ge. gamsq) go to 220
               assign 210 to igo
c              fix-h..
               go to 70
  210          continue
               sd2=sd2/gam**2
               sh21=sh21*gam
               sh22=sh22*gam
          go to 200
  220 continue
          if(sflag)250,230,240
  230     continue
               sparam(3)=sh21
               sparam(4)=sh12
               go to 260
  240     continue
               sparam(2)=sh11
               sparam(5)=sh22
               go to 260
  250     continue
               sparam(2)=sh11
               sparam(3)=sh21
               sparam(4)=sh12
               sparam(5)=sh22
  260 continue
          sparam(1)=sflag
          return
      end
