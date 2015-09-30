!.Structured grid module
module modMSH

  use modGLB
  implicit real(dbl_kind)(a-h,o-z), integer(i-n)

!.variables
  public
  integer                                         :: IN, JN     ! # of grid in x dir
  real(dbl_kind)                                   :: alpha, gamma  
  real(dbl_kind)                                   :: dx, dy    
  character(len=30)                               :: cfGrd      
  
!.basic varibles of the geometry
  real(dbl_kind), dimension(:,:), allocatable :: xnd, ynd, znd   ! position of grid node 
  real(dbl_kind), dimension(:,:), allocatable :: xcl, ycl, zcl   ! position of cell
  real(dbl_kind), dimension(:,:), allocatable :: xsv, ysv, zsv   ! position of vertical side  
  real(dbl_kind), dimension(:,:), allocatable :: xsh, ysh, zsh   ! position of horizontal side

!.function
  public  :: sbmRead   ! read file 
  private :: sbmInti   ! initalization  
  private :: sbmAllc   ! allocate module variables 

  contains
  
  subroutine sbmRead(IDFile, cFile)
  
    use modGLB
    implicit real(dbl_kind)(a-h,o-z), integer(i-n) 

    integer,           intent(in) :: IDFile
    character(len=25), intent(in) :: cFile
    character(len=15)              :: KEYWORD
 
    logical            :: exFile   
    integer            :: IDGrd = 8

    read(IDFile, '(A2, A13, A15)', err=910, end=950) BLANK, KEYWORD, cfGrd
	write(21, *) KEYWORD, cfGrd
    if(KEYWORD/='GridFile') goto 930

    read(IDFile, *, err=910, end=950) KEYWORD, iN
    write(21, *) KEYWORD, iN
    if(KEYWORD/='iN') goto 930
    
    read(IDFile, *, err=910, end=950) KEYWORD, jN
    write(21, *) KEYWORD, jN
    if(KEYWORD/='jN') goto 930

    read(IDFile, *, err=910, end=950) KEYWORD, dx
    write(21, *) KEYWORD, dx
    if(KEYWORD/='dx') goto 930
    
    read(IDFile, *, err=910, end=950) KEYWORD, dy
    write(21, *) KEYWORD, dy
    if(KEYWORD/='dy') goto 930    

    call sbmAllc
    
!    alpha = 0.006d0;   gamma = 2.0d0    
    alpha = 1.0d0;     gamma = 1.0d0    
    
    select case(cfGrd)
!     case('Rec20')
!      do j=-1, jN+1
!        do i=-1, iN+1
!          xnd(i,j) = (dble(i)-0.5d0)*dx
!          ynd(i,j) = (dble(j)-0.5d0)*dy 
!          znd(i,j) = - 1.0d0
!        end do
!      end do

!     case('WetDry')
!      beta1 = 0.001d0;   beta2 = 0.01d0
!      do j=-1, jN+1
!        do i=-1, iN+1
!          xnd(i,j) = (dble(i)-0.5d0)*dx
!          ynd(i,j) = (dble(j)-0.5d0)*dy 
!          if(xnd(i,j) .le. 100.0d0) then
!            znd(i,j) = - beta1*xnd(i,j)
!          else if(xnd(i,j) .gt. 100.0d0 .and. xnd(i,j) .le. 200.0d0) then
!            znd(i,j) = - 0.1d0 - beta2*(xnd(i,j)-100.0d0)
!          else if(xnd(i,j) .gt. 200.0d0) then
!            znd(i,j) = - 1.1d0 - beta1*(xnd(i,j)-200.0d0)
!          end if
!          znd(i,j) = znd(i,j) + 1.40d0
!          znd(i,j) = 1.4d0-0.005d0*xnd(i,j)
!        end do
!      end do
      
     case('shear')   ! dimensional, dx=0.01m,dy=0.01m
      do j=-1, jN+1
        do i=-1, iN+1
          xnd(i,j) = dble(i)*dx
          ynd(i,j) = dble(j)*dy 
          znd(i,j) = -0.066d0   !case VI, here let znd be positive, different from equation, usually znd is negative
        end do 
      end do     

!     case default
!      open(IDGrd, file=cfGrd, status='Old')
!      read(IDGrd, *, err=101, end=102) iN, jN
!      iN = iN - 2;   jN = jN - 2
!      read(IDGrd, *, err=101, end=102) iGE, jGN, iGW, jGS
!      do j=-1, jN+1
!        do i=-1, iN+1
!          read(IDGrd,*) xnd(i,j), ynd(i,j), znd(i,j)
!          znd(i,j) = -znd(i,j)
!        end do
!      end do
!      close(IDGrd)      
    
    end select
    
    call sbmInit
    
    return
    
101 write(21,*) 'FATAL: error reading data from grid file ', cFile
    stop
102 write(21,*) 'FATAL: unexpected end of grid file ', cFile
    stop    
    
910 write(21,*) 'FATAL: error reading control file,', cFile
    stop
930 write(21,*) 'FATAL: KEYWORD does not defined!', KEYWORD
    stop
950 write(21,*) 'FATAL: unexpect end of control file,', cFile 
    stop    
    
  end subroutine sbmRead  
  
  
!.Initialization and grid step
  subroutine sbmInit

    implicit real(dbl_kind)(a-h,o-z), integer(i-n)

  !.calculate cell center position 
    do i=0, iN+1
      do j=0, jN+1
        xcl(i,j) = 0.5d0 * (xnd(i,j) + xnd(i-1,j))
        ycl(i,j) = 0.5d0 * (ynd(i,j) + ynd(i,j-1))
        zcl(i,j) = znd(i,j)  ! since znd is const for shear flow
      end do
    end do
  !.linear calculation of side position
    do i=-1, iN+1
      do j=0, jN+1
        xsv(i,j) = xnd(i,j)
        ysv(i,j) = 0.5d0 * (ynd(i,j) + ynd(i,j-1))
        zsv(i,j) = znd(i,j)
      end do
    end do
    do i=0, iN+1
      do j=-1, jN+1
        xsh(i,j) = 0.5d0 * (xnd(i,j) + xnd(i-1,j))
        ysh(i,j) = ynd(i,j)
        zsh(i,j) = znd(i,j)
      end do
    end do

    open (100, file='ShearGrd.dat')
    write(100, *) 'Title="Grid"'
    write(100, *) 'Variables = "x", "y", "z"' 
    write(100, *) 'Zone T="GridPlot"'
    write(100, '(A3, I5, A6, I5, A11)') ' I=', iN+2, ',   J=', jN+2, ',   f=point'
    do j=0, jN+1
      do i=0, iN+1
        write(100, '(8E30.20)') xcl(i,j), ycl(i,j), zcl(i,j) 
      end do
    end do
    
    return

  end subroutine sbmInit



!//allocate memory for variables
  subroutine sbmAllc
  !.allocation
    iStat = 0
  !.node information
    if(iStat==0) allocate(xnd(-1:IN+1, -1:jN+1), stat=iStat)
    if(iStat==0) allocate(ynd(-1:IN+1, -1:jN+1), stat=iStat)
    if(iStat==0) allocate(znd(-1:IN+1, -1:jN+1), stat=iStat)

  !.vertical side information
    if(iStat==0) allocate(xsv(-1:IN+1, 0:jN+1), stat=iStat)
    if(iStat==0) allocate(ysv(-1:IN+1, 0:jN+1), stat=iStat)
    if(iStat==0) allocate(zsv(-1:IN+1, 0:jN+1), stat=iStat)
  !.horizontal side information
    if(iStat==0) allocate(xsh(0:IN+1, -1:jN+1), stat=iStat)
    if(iStat==0) allocate(ysh(0:IN+1, -1:jN+1), stat=iStat)
    if(iStat==0) allocate(zsh(0:IN+1, -1:jN+1), stat=iStat)
  !.cell information
    if(iStat==0) allocate(xcl(0:IN+1, 0:jN+1), stat=iStat)
    if(iStat==0) allocate(ycl(0:IN+1, 0:jN+1), stat=iStat)
    if(iStat==0) allocate(zcl(0:IN+1, 0:jN+1), stat=iStat)

  !.detect allocation problem 
    if(iStat/=0) then
      write(21,*) 'FATAL: allocation problem "subGrdAlc"!'
      write(21,*) 'Program terminates abnormally!'
      stop   
    end if

    xnd = 0.0d0;   ynd = 0.0d0;   znd = 0.0d0
    xsv = 0.0d0;   ysv = 0.0d0;   zsv = 0.0d0
    xsh = 0.0d0;   ysh = 0.0d0;   zsh = 0.0d0
    xcl = 0.0d0;   ycl = 0.0d0;   zcl = 0.0d0    

    return

  end subroutine sbmAllc

end module modMSH
