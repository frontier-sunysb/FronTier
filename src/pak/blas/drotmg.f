      subroutine drotmg (dd1,dd2,dx1,dy1,dparam)
c
c     construct the modified givens transformation matrix h which zeros
c     the second component of the 2-vector  (dsqrt(dd1)*dx1,dsqrt(dd2)*
c     dy2)**t.
c     with dparam(1)=dflag, h has one of the following forms..
c
c     dflag=-1.d0     dflag=0.d0        dflag=1.d0     dflag=-2.d0
c
c       (dh11  dh12)    (1.d0  dh12)    (dh11  1.d0)    (1.d0  0.d0)
c     h=(          )    (          )    (          )    (          )
c       (dh21  dh22),   (dh21  1.d0),   (-1.d0 dh22),   (0.d0  1.d0).
c     locations 2-4 of dparam contain dh11, dh21, dh12, and dh22
c     respectively. (values of 1.d0, -1.d0, or 0.d0 implied by the
c     value of dparam(1) are not stored in dparam.)
c
c     the values of gamsq and rgamsq set in the data statement may be
c     inexact.  this is ok as they are only used for testing the size
c     of dd1 and dd2.  all actual scaling of data is done using gam.
c
      double precision gam,one,rgamsq,dd2,dh11,dh21,dparam,dp2,
     1 dq2,du,dy1,zero,gamsq,dd1,dflag,dh12,dh22,dp1,dq1,
     2 dtemp,dx1,two
      dimension dparam(5)
c
      data zero,one,two /0.d0,1.d0,2.d0/
      data gam,gamsq,rgamsq/4096.d0,16777216.d0,5.9604645d-8/
      if(.not. dd1 .lt. zero) go to 10
c       go zero-h-d-and-dx1..
          go to 60
   10 continue
c     case-dd1-nonnegative
      dp2=dd2*dy1
      if(.not. dp2 .eq. zero) go to 20
          dflag=-two
          go to 260
c     regular-case..
   20 continue
      dp1=dd1*dx1
      dq2=dp2*dy1
      dq1=dp1*dx1
c
      if(.not. dabs(dq1) .gt. dabs(dq2)) go to 40
          dh21=-dy1/dx1
          dh12=dp2/dp1
c
          du=one-dh12*dh21
c
          if(.not. du .le. zero) go to 30
c         go zero-h-d-and-dx1..
               go to 60
   30     continue
               dflag=zero
               dd1=dd1/du
               dd2=dd2/du
               dx1=dx1*du
c         go scale-check..
               go to 100
   40 continue
          if(.not. dq2 .lt. zero) go to 50
c         go zero-h-d-and-dx1..
               go to 60
   50     continue
               dflag=one
               dh11=dp1/dp2
               dh22=dx1/dy1
               du=one+dh11*dh22
               dtemp=dd2/du
               dd2=dd1/du
               dd1=dtemp
               dx1=dy1*du
c         go scale-check
               go to 100
c     procedure..zero-h-d-and-dx1..
   60 continue
          dflag=-one
          dh11=zero
          dh12=zero
          dh21=zero
          dh22=zero
c
          dd1=zero
          dd2=zero
          dx1=zero
c         return..
          go to 220
c     procedure..fix-h..
   70 continue
      if(.not. dflag .ge. zero) go to 90
c
          if(.not. dflag .eq. zero) go to 80
          dh11=one
          dh22=one
          dflag=-one
          go to 90
   80     continue
          dh21=-one
          dh12=one
          dflag=-one
   90 continue
      go to igo,(120,150,180,210)
c     procedure..scale-check
  100 continue
  110     continue
          if(.not. dd1 .le. rgamsq) go to 130
               if(dd1 .eq. zero) go to 160
               assign 120 to igo
c              fix-h..
               go to 70
  120          continue
               dd1=dd1*gam**2
               dx1=dx1/gam
               dh11=dh11/gam
               dh12=dh12/gam
          go to 110
  130 continue
  140     continue
          if(.not. dd1 .ge. gamsq) go to 160
               assign 150 to igo
c              fix-h..
               go to 70
  150          continue
               dd1=dd1/gam**2
               dx1=dx1*gam
               dh11=dh11*gam
               dh12=dh12*gam
          go to 140
  160 continue
  170     continue
          if(.not. dabs(dd2) .le. rgamsq) go to 190
               if(dd2 .eq. zero) go to 220
               assign 180 to igo
c              fix-h..
               go to 70
  180          continue
               dd2=dd2*gam**2
               dh21=dh21/gam
               dh22=dh22/gam
          go to 170
  190 continue
  200     continue
          if(.not. dabs(dd2) .ge. gamsq) go to 220
               assign 210 to igo
c              fix-h..
               go to 70
  210          continue
               dd2=dd2/gam**2
               dh21=dh21*gam
               dh22=dh22*gam
          go to 200
  220 continue
          if(dflag)250,230,240
  230     continue
               dparam(3)=dh21
               dparam(4)=dh12
               go to 260
  240     continue
               dparam(2)=dh11
               dparam(5)=dh22
               go to 260
  250     continue
               dparam(2)=dh11
               dparam(3)=dh21
               dparam(4)=dh12
               dparam(5)=dh22
  260 continue
          dparam(1)=dflag
          return
      end
