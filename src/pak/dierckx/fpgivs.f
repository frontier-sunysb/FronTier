      subroutine fpgivs(piv,ww,cos,sin)
      implicit none
c  subroutine fpgivs calculates the parameters of a givens
c  transformation .
c  ..
c  ..scalar arguments..
      double precision piv,ww,cos,sin
c  ..local scalars..
      double precision dd,one,store
c  ..function references..
      double precision abs,dsqrt
c  ..
      one = 0.1d+01
      store = abs(piv)
      if(store.ge.ww) dd = store*dsqrt(one+(ww/piv)**2)
      if(store.lt.ww) dd = ww*dsqrt(one+(piv/ww)**2)
      cos = ww/dd
      sin = piv/dd
      ww = dd
      return
      end
