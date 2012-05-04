      subroutine fprota(cos,sin,a,b)
      implicit none
c  subroutine fprota applies a givens rotation to a and b.
c  ..
c  ..scalar arguments..
      double precision cos,sin,a,b
c ..local scalars..
      double precision stor1,stor2
c  ..
      stor1 = a
      stor2 = b
      b = cos*stor2+sin*stor1
      a = cos*stor1-sin*stor2
      return
      end
