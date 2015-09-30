!//define global parameters
module modGLB
!.variables
  integer, parameter        :: dbl_kind = 8             !  type of double parmeter
  integer, parameter        :: IDB = 0                  !  debug mode
  real(dbl_kind), parameter  :: third = 1.0d0/3.0d0
  real(dbl_kind), parameter  :: grv = 9.81d0             !  gravity accerleration
  real(dbl_kind), parameter  :: rho = 1025.0d0           !  water density
  real(dbl_kind), parameter  :: hTol = 1.0E-6            !  tolerant water depth

  real(dbl_kind)  :: pi                                   !

  contains
    
  subroutine constant(pii)
    
    implicit real(dbl_kind)(a-h,o-z), integer(i-n)
    
    pi = dacos(-1.0d0)
  
    return
    
  end subroutine
  


!.newton's method for calculate wavenumber 
  subroutine newton(aDep, aPrd, aWk)
  
    implicit real(dbl_kind)(a-h,o-z), integer(i-n)  
  
    omg = 2.0d0*pi/aPrd
    wk0 = omg/dsqrt(grv*aDep)
    do it = 1, 1000
      wkh  = wk0*aDep
      fnct = omg*omg - grv*wk0*dtanh(wkh)
      dfdk = - grv*(dtanh(wkh)+wkh/dcosh(wkh)**2.0d0)
      aWk  = wk0 - fnct/dfdk
      if(dabs(aWk-wk0) .le. 1.0E-20) exit 
      wk0 = aWk
    end do 
       
    return

  end subroutine  
    
  
  
end module