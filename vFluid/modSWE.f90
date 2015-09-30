module modSWE
  
  use modGLB;   use modMSH  
  implicit real(dbl_kind)(a-h,o-z), integer(i-n)

!.variables
  real(dbl_kind)     :: cfbm   !.bottom friction coefficient  
  real(dbl_kind)     :: cf     !.friction coefficient
  real(dbl_kind)     :: l      !.velocity profile parameter
  character(len=3) :: cAdv
!.sparse matrix solution choices for itpack
  integer           :: iSlv   !.sparse matrix solver(jcg,jsi,sor...)     
  integer           :: mItr   !.maximum iteration, 1000
  integer           :: iRmv   !.remove, 0
  real(dbl_kind)     :: theta  !.theta
  real(dbl_kind)     :: zeta   !.stopping criterion, 1.0e-6
  real(dbl_kind)     :: tolr   !.tolerance, 1.0e-12  
!.boundary condition  
  character(len=3) :: cSweBE, cSweBS, cSweBW, cSweBN   ! east/north boundary type (InnFlw/OutFlw/SlpWll)
  real(dbl_kind)     :: etaBE, qxcBE, qycBE             ! prescribed values at east boundary
  real(dbl_kind)     :: etaBN, qxcBN, qycBN             ! prescribed values at north boundary
  real(dbl_kind)     :: etaBW, qxcBW, qycBW             ! prescribed values at west boundary
  real(dbl_kind)     :: etaBS, qxcBS, qycBS             ! prescribed values at south boundary 
!.primitive varibles
  real(dbl_kind), public, dimension(:,:), allocatable :: eta, dpc   ! water level and depth
  real(dbl_kind), public, dimension(:,:), allocatable :: usd, vsd   ! velocity at side
  real(dbl_kind), public, dimension(:,:), allocatable :: ucl, vcl   ! velocity at cell center
  real(dbl_kind), public, dimension(:,:), allocatable :: dsv, dsh     ! depth at cell side  
  real(dbl_kind), dimension(:,:), allocatable :: etan, usdn, vsdn    ! new step
  real(dbl_kind), dimension(:,:), allocatable :: und, vnd      !velocity at node
!.update velocity due to advection 
  real(dbl_kind), dimension(:,:), allocatable :: ustr, vstr   ! advection term
!.velocity gradient at cell center 
  real(dbl_kind) ::  dudx, dudy     ! dudx, dudy at cell center
  real(dbl_kind) ::  dvdx, dvdy     ! dvdx, dvdy at cell center
!.intermediate variables 
  real(dbl_kind), dimension(:,:), allocatable :: xBfs, yBfs   ! bottom friction coeficient
  real(dbl_kind), dimension(:,:), allocatable :: xGht, yGht   ! sum of explicit terms Ghat
!.sparse matrix related
  integer,       dimension(:,:,:), allocatable :: iMap       !.map ix,jy,ij -> solution vector 
  real(dbl_kind), dimension(:,:,:), allocatable :: spm        !.matrix coefficient formed in (ix,jy,ij)
  real(dbl_kind), dimension(:,:),   allocatable :: arsv       !.inverse of A at vertical side
  real(dbl_kind), dimension(:,:),   allocatable :: arsh       !.inverse of A at horizontal side 
  real(dbl_kind), dimension(:,:),   allocatable :: rhs        !.R.H.S of the continuity equations
  real(dbl_kind), dimension(:,:), allocatable   :: Umean, Vmean, Um, Vm

!.function
  public  :: sbcREAD      ! read 
  public  :: sbcINIT      ! initialization
  public  :: sbcSOLV      ! solve precedure  
  private :: sbcALLC      ! allocate module variables
  private :: sbcGHT       ! explicit terms 
  private :: sbcSPM       ! assemble the sparse matrix M
  private :: sbcETA       ! update the water elevation
  public :: sbcBCUPD     ! numerical boundary condition

contains

!.initialization
  subroutine sbcREAD(IDFile, cFile)
  
    use modGLB
    implicit real(dbl_kind)(a-h,o-z), integer(i-n)  

    integer,           intent(in)  :: IDFile
    character(len=15), intent(in)  :: cFile
    character(len=15)              :: KEYWORD

    write(21,*) ' call sbcRead... '

  !.variable allocation according to mesh dimension
    call sbcALLC

  !.parameters reading 
    read(IDFile, *, err=910, end=950) KEYWORD, cfbm
    write(21, *) KEYWORD, cfbm
    if(KEYWORD/='BtmFriction') goto 930
    
    read(IDFile, *, err=910, end=950) KEYWORD, cAdv
    write(21, *) KEYWORD, cAdv
    if(KEYWORD/='Advection') goto 930
    
    read(IDFile, *, err=910, end=950) KEYWORD, cSweBE
    write(21, *) KEYWORD, cSweBE
    if(KEYWORD/='cSweBE') goto 930
    if(cSweBE=='Dis') then
      read(IDFile, *, err=910, end=950) KEYWORD, qxcBE
      write(21, *) KEYWORD, qxcBE
      if(KEYWORD/='qxcBE') goto 930
      read(IDFile, *, err=910, end=950) KEYWORD, qycBE
      write(21, *) KEYWORD, qycBE
      if(KEYWORD/='qycBE') goto 930		  	
    else if(cSweBE=='Eta') then
      read(IDFile, *, err=910, end=950) KEYWORD, etaBE
      write(21, *) KEYWORD, etaBE
      if(KEYWORD/='etaBE') goto 930
    end if
  
    read(IDFile, *, err=910, end=950) KEYWORD, cSweBW
    write(21,*) KEYWORD, cSweBW
    if(KEYWORD/='cSweBW') goto 930
    if(cSweBW=='Dis') then
      read(IDFile, *, err=910, end=950) KEYWORD, qxcBW
      write(21,*) KEYWORD, qxcBW
      if(KEYWORD/='qxcBW') goto 930
      read(IDFile, *, err=910, end=950) KEYWORD, qycBW
      write(21,*) KEYWORD, qycBW
      if(KEYWORD/='qycBW') goto 930
    else if(cSweBW=='Eta') then
      read(IDFile, *, err=910, end=950) KEYWORD, etaBW
      write(21,*) KEYWORD, etaBW
      if(KEYWORD/='etaBW') goto 930
    end if
    
    read(IDFile, *, err=910, end=950) KEYWORD, cSweBS
    write(21,*) KEYWORD, cSweBS
    if(KEYWORD/='cSweBS') goto 930
    if(cSweBS=='Dis') then
      read(IDFile, *, err=910, end=950) KEYWORD, qxcBS
      write(21,*) KEYWORD, qxcBS
      if(KEYWORD/='qxcBS') goto 930
       read(IDFile, *, err=910, end=950) KEYWORD, qycBS
       write(21,*) KEYWORD, qycBS
      if(KEYWORD/='qycBS') goto 930		  	
    else if(cSweBS=='Eta') then
      read(IDFile, *, err=910, end=950) KEYWORD, etaBS
      write(21,*) KEYWORD, etaBS
      if(KEYWORD/='etaBS') goto 930
    end if    

    read(IDFile, *, err=910, end=950) KEYWORD, cSweBN
    write(21,*) KEYWORD, cSweBN
    if(KEYWORD/='cSweBN') goto 930
    if(cSweBN=='Dis') then
      read(IDFile, *, err=910, end=950) KEYWORD, qxcBN
      write(21,*) KEYWORD, qxcBN
      if(KEYWORD/='qxcBN') goto 930
      read(IDFile, *, err=910, end=950) KEYWORD, qycBN
      write(21,*) KEYWORD, qycBN
      if(KEYWORD/='qycBN') goto 930
    else if(cSweBN=='Eta') then
      read(IDFile, *, err=910, end=950) KEYWORD, etaBN
      write(21,*) KEYWORD, etaBN
      if(KEYWORD/='etaBN') goto 930
    end if

    read(IDFile, *, err=910, end=950) KEYWORD, theta
    write(21,*) KEYWORD, theta        
    if(KEYWORD/='Theta') goto 930
       
    read(IDFile, *, err=910, end=950) KEYWORD, iSlv
    write(21,*) KEYWORD, iSlv    
    if(KEYWORD/='SpmSolver') goto 930 
    
    read(IDFile, *, err=910, end=950) KEYWORD, mItr
    write(21,*) KEYWORD, mItr        

    if(KEYWORD/='MaxIter') goto 930
    
    read(IDFile, *, err=910, end=950) KEYWORD, iRmv
    write(21,*) KEYWORD, iRmv    
    if(KEYWORD/='Remove') goto 930   
    
    read(IDFile, *, err=910, end=950) KEYWORD, zeta
    write(21,*) KEYWORD, zeta        
    if(KEYWORD/='Zeta') goto 930
    
    read(IDFile, *, err=910, end=950) KEYWORD, tolr
    write(21,*) KEYWORD, tolr        
    if(KEYWORD/='Tolerance') goto 930  

    write(21,*) ' end call sbcRead... '

    return

910 write(21,*) 'FATAL: error reading control file,', cFile
    stop
930 write(21,*) 'FATAL: KEYWORD does not defined!', KEYWORD
    stop
950 write(21,*) 'FATAL: unexpect end of control file,', cFile 
    stop

  end subroutine sbcREAD
  


!.initialization
  subroutine sbcINIT(iHot, IDHot, minStp, tme)
  
    implicit real(dbl_kind)(a-h,o-z), integer(i-n)
!.hot start   
!    real(dbl_kind), dimension(:,:), allocatable :: eta00
!    real(dbl_kind), dimension(:,:), allocatable :: usd00
!    real(dbl_kind), dimension(:,:), allocatable :: vsd00

    write(21,*) ' call sbcInit... '

  !.Initialization or cold start
    minStp = 1;   tme = 0.0d0
  !.primitive variables   
    eta = 0.0d0;   dpc = 0.0d0
    usd = 0.0d0;   vsd = 0.0d0
    ucl = 0.0d0;   vcl = 0.0d0 
    dsv = 0.0d0;   dsh = 0.0d0  
    ustr = 0.0d0;  vstr = 0.0d0
  !.force terms  
    xBfs = 0.0d0;  xGht = 0.0d0 
    yBfs = 0.0d0;  yGht = 0.0d0 
    Umean = 0.0d0; Vmean = 0.0d0
    Um = 0.0d0; Vm = 0.0d0
  !.matrix related terms  
    arsv = 0.0d0;  arsh = 0.0d0 
    imap = 0;      spm = 0.0d0;   rhs = 0.0d0  
!.hot start  
!    if(iHot==1) then
!       read(IDHot, err=910) minStp, tme, IN0, JN0
!       if(IN0/=iN .or. JN0/=jN) then
!         write(21,*) " Warning: dimension miss match in restart file" 
!         istat = 0
!         if(istat==0) allocate(eta00(0:iN0+1,0:jN0+1), stat=istat)
!         if(istat==0) allocate(usd00(0:iN0+1,0:jN0+1), stat=istat)
!         if(istat==0) allocate(vsd00(0:iN0+1,0:jN0+1), stat=istat)
!       end if
!	   read (IDHot, err=910) ((eta00(i,j), i=0, iN0+1), j=0, jN0+1),   &
!                     &        ((usd00(i,j), i=0, iN0+1), j=0, jN0+1),   &
!                     &        ((vsd00(i,j), i=0, iN0+1), j=0, jN0+1)    
!       do i=0, iN+1
!         do j=0, jN+1
!           eta(i,j) = eta00(i,5)
!           usd(i,j) = usd00(i,5)
!           vsd(i,j) = vsd00(i,5)
!           usd(i,j) = 0.0d0;   vsd(i,j) = 0.0d0
!         end do
!       end do
!       if(allocated(eta00)) deallocate(eta00)
!       if(allocated(usd00)) deallocate(usd00)
!       if(allocated(vsd00)) deallocate(vsd00)       
!    end if

  !.The previous step water elevation
    etan = eta;   usdn = usd;   vsdn = vsd

  !.Intialize boundary values
    call sbcBCUPD(tme+dtme,dtme)
  !.interpolation
    call sbcAuxVar
  !.Connectivity between cell and equations
    iMap = 0
    do jy=1, jN
      do ix=1, iN
        iMap(ix,jy,0) = iN*(jy-1) + ix       ! ix,jy
        iMap(ix,jy,1) = iN*(jy-1) + ix + 1   ! ix+1,jy
        iMap(ix,jy,2) = iN*jy + ix           ! ix,jy+1
        iMap(ix,jy,3) = iN*(jy-1) + ix - 1   ! ix-1,jy
        iMap(ix,jy,4) = iN*(jy-2) + ix       ! ix,jy-1
      end do
    end do

    write(21,*) ' end call sbcInit... '

    return
    
910 write(21,*) 'FATAL: error reading restart file,'
    stop
950 write(21,*) 'FATAL: unexpect end of restart file,' 
    stop

  end subroutine sbcINIT
  
  
  
!.allocate variables 
  subroutine sbcALLC
  
    use modGLB
    implicit real(dbl_kind)(a-h,o-z), integer(i-n)

    write(21,*) ' call sbcAllc... '

    iStat = 0  
  !.allocation for primitive variables
    if(iStat==0) allocate(eta(0:iN+1, 0:jN+1), stat=iStat)
    if(iStat==0) allocate(dpc(1:iN, 1:jN), stat=iStat)
    if(iStat==0) allocate(ucl(0:iN+1, 0:jN+1), stat=iStat)
    if(iStat==0) allocate(vcl(0:iN+1, 0:jN+1), stat=iStat)
    if(iStat==0) allocate(usd(-1:iN+1, 0:jN+1), stat=iStat) !-1 and 0 are ghost point u(ix,1) is in the domain
    if(iStat==0) allocate(vsd(0:iN+1, -1:jN+1), stat=iStat) !-1 and jN+1 are ghost point, boundary is 0 and jN
  !.next step  
    if(iStat==0) allocate(etan(0:iN+1, 0:jN+1), stat=iStat)
    if(iStat==0) allocate(ustr(0:iN, 1:jN+1), stat=iStat)
    if(iStat==0) allocate(usdn(1:iN, 1:jN+1), stat=iStat)
    if(iStat==0) allocate(vstr(1:iN+1, 0:jN), stat=iStat)
    if(iStat==0) allocate(vsdn(1:iN+1, 1:jN), stat=iStat)
    if(iStat==0) allocate(Umean(0:iN,1:jN), stat=iStat)
    if(iStat==0) allocate(Vmean(0:iN,1:jN), stat=iStat)
    if(iStat==0) allocate(Um(0:iN,1:jN), stat=iStat)
    if(iStat==0) allocate(Vm(0:iN,1:jN), stat=iStat)
  !.allocation for side information
    if(iStat==0) allocate(und(0:iN, 0:jN), stat=iStat)
    if(iStat==0) allocate(vnd(0:iN, 0:jN), stat=iStat)
  !.allocation for side information
    if(iStat==0) allocate(dsv(0:iN, 1:jN+1), stat=iStat)
    if(iStat==0) allocate(dsh(1:iN+1, 0:jN), stat=iStat)
  !.allocatation for source
    if(iStat==0) allocate(xGht(0:iN, 1:jN+1), stat=iStat)
    if(iStat==0) allocate(xBfs(0:iN, 1:jN+1), stat=iStat)
    if(iStat==0) allocate(yGht(1:iN+1, 0:jN), stat=iStat)
    if(iStat==0) allocate(yBfs(1:iN+1, 0:jN), stat=iStat)
  !.allocatation for sparse matrix related
    if(iStat==0) allocate(arsv(0:iN, 1:jN+1), stat=iStat)
    if(iStat==0) allocate(arsh(1:iN+1, 0:jN), stat=iStat)
    if(iStat==0) allocate(rhs (1:iN, 1:jN),     stat=iStat)
    if(iStat==0) allocate(iMap(1:iN, 1:jN, 0:4),stat=iStat)
    if(iStat==0) allocate(spm (1:iN, 1:jN, 0:4),stat=iStat)
   
    if(iStat/=0) then
      write(21,*) ' FATAL: allocation problem "sbcVarAlc"!'
      write(21,*) ' Program terminates abnormally!'
      stop   
    end if   

    write(21,*) ' end call sbcAllc... '

    return

  end subroutine sbcALLC



  subroutine sbcSOLV(tme, dtme)

    implicit real(dbl_kind)(a-h,o-z), integer(i-n) 
  
  !.calculate the explicit terms
    call sbcGHT(tme, dtme)
  !.update the water elevation  
    call sbcETA(tme, dtme)
  !.update the velocity
    do jy=1, jN  ! try update to jN
      do ix=1, iN-1
        ip =  ix +1
        eta_x  = (eta(ip,jy)-eta(ix,jy))/dx
        usdn(ix,jy) = xGht(ix,jy) - dtme*grv*dsv(ix,jy)*eta_x*arsv(ix,jy)

      end do
    end do

    do ix=2, iN  !try update to iN
      do jy=1, jN-1
        jp =  jy + 1  ! backward for eta, too
        eta_y = (eta(ix,jp)-eta(ix,jy))/dy

        vsdn(ix,jy) = yGht(ix,jy) - dtme*grv*dsh(ix,jy)*eta_y*arsh(ix,jy)

      ! if(vsdn(ix,jy) .LE. 1.0E-16) vsdn(ix,jy) = 0.0d0

      end do
    end do

    usd = usdn;   vsd = vsdn
  !.update boundary values
    call sbcBCUPD(tme+dtme,dtme)
  !.interpolation
    call sbcAuxVar
    
    return

  end subroutine sbcSOLV



  subroutine sbcAuxVar
  
    implicit real(dbl_kind)(a-h,o-z), integer(i-n)
  !.calculate depth of water at cell center
    do ix=1, iN
      do jy=1, jN
        dpc(ix,jy) = eta(ix,jy) - zcl(ix,jy)  !h=d+eta
      end do
    end do
                
  !.calculate dsv
    do ix=0, iN
      do jy=1, jN

       ! if(ix==-1) etas = eta(ix+1,jy)
       ! if(ix==iN+1) etas = eta(iN,jy)
       ! if(ix/=-1 .and. ix/=iN+1) etas = 0.5d0*(eta(ix,jy) + eta(ix+1,jy))
        etas = 0.5d0*(eta(ix,jy) + eta(ix+1,jy))
        dsv(ix,jy) = etas - zsv(ix,jy)
      end do
    end do

    do ix=1, iN
      do jy=0, jN

       ! if(jy==-1) etas = eta(ix,jy+1)
       ! if(jy==jN+1) etas = eta(ix,jy-1)
       ! if(jy/=-1 .and. jy/=jN+1) etas = 0.5d0*(eta(ix,jy) + eta(ix,jy+1))
        etas = 0.5d0*(eta(ix,jy) + eta(ix,jy+1))
        dsh(ix,jy) = etas - zsh(ix,jy)
      end do
    end do
  !.calculate velocity at mesh node
    und = 0.0d0;   vnd = 0.0d0
    do jy=0, jN
      do ix=0, iN
        und(ix,jy) = 0.5d0 * (usd(ix,jy)+usd(ix,jy+1))
        vnd(ix,jy) = 0.5d0 * (vsd(ix,jy)+vsd(ix+1,jy))
      end do
!    und(iN+1,jy) = und(iN,jy);   vnd(iN+1,jy) = vnd(iN,jy)
!    und(-1,jy) = und(0,jy);   vnd(-1,jy) = 0
    end do
    do ix=0, iN
!    und(ix,jN+1) = und(ix,jN);   vnd(ix,jN+1) = 0
    und(ix,0) = und(ix,1);   vnd(ix,0) = vnd(ix,1)
    end do

  !.calculate velocity at cell center
    ucl = 0.0d0;   vcl = 0.0d0
    do ix=0, iN+1
      do jy=0, jN+1
        ucl(ix,jy) = 0.5d0 * (usd(ix-1,jy)+usd(ix,jy))
        vcl(ix,jy) = 0.5d0 * (vsd(ix,jy-1)+vsd(ix,jy))
      end do

    end do

    return
  
  end subroutine 


  subroutine sbcGHT(tme, dtme)

    implicit real(dbl_kind)(a-h,o-z), integer(i-n)
 
  !.Advection Terms
!    call sbcCDS(dtme)
    call sbcUPW(dtme)   !calculate ustr,vstr
  !.Bottom Friction Coefficient,at cell side
    xBfs = 0.0d0;   yBfs = 0.0d0
    do jy=1, jN
      do ix=0, iN
    
         if(jy.LE.third*jN) then
          cf = cfbm + 28.5d0*dsv(ix,jy) !bottom friction + vegetated drag
         else
          cf = cfbm
         end if
         us = 0.5d0*(ucl(ix,jy) + ucl(ix+1,jy));  vs = vsd(ix,jy)
         absVel = dsqrt(us**2.0d0 + vs**2.0d0)
         xBfs(ix,jy) = 0.5d0*cf*absVel

      end do
    end do

    do ix=1, iN
      do jy=0, jN
 
         if(jy.LE.third*jN) then
          cf = cfbm + 28.5d0*dsh(ix,jy) !bottom friction + vegetated drag
         else
          cf = cfbm
         end if
         us = usd(ix,jy);  vs = 0.5d0*(vcl(ix,jy) + vcl(ix,jy+1))
         absVel = dsqrt(us**2.0d0 + vs**2.0d0)
         yBfs(ix,jy) = 0.5d0*cf*absVel


      end do
    end do
    
  !.Combine the Explicit Terms in Ghat
    xGht = 0.0d0;   yGht = 0.0d0
    do jy=1, jN
      do ix=0, iN

        arsv(ix,jy) = 1.0d0 / (dsv(ix,jy) + dtme*xBfs(ix,jy))
        xGht(ix,jy) = arsv(ix,jy)*ustr(ix,jy)*dsv(ix,jy)

      end do
    end do

    do ix=1, iN
      do jy=0, jN

        arsh(ix,jy) = 1.0d0 / (dsh(ix,jy) + dtme*yBfs(ix,jy))
        yGht(ix,jy) = arsh(ix,jy)*vstr(ix,jy)*dsh(ix,jy)

!        if(dabs(yGht(ix,jy))>1.0E-10) then
!          write(6,*) ix,jy
!          write(6,*) vstr(ix,jy), arsh(ix,jy), dsh(ix,jy), yRfs(ix,jy), yDfs(ix,jy), ybfs(ix,jy)
!          pause
!        end if       
        
      end do
    end do

    return

  end subroutine sbcGHT



  subroutine sbcUPW(dtme)

    implicit real(dbl_kind)(a-h,o-z), integer(i-n)
      
    ustr = 0.0d0;   vstr = 0.0d0
    do jy=1, jN
      do ix=0, iN
      
      !.upwinding difference
        us = usd(ix,jy)
        iu = ix - int(max(0.,dsign(1.1d0,us)))   !us>0,so iu=ix-1, backward
        ic = iu + 1  !ix
        rdx1 = 1.0d0/(xsv(ic,jy)-xsv(iu,jy))
        phi_x = 1.0d0/(alpha*gamma*xsv(ix,jy)**(gamma-1.0d0))
        dudx = phi_x*(usd(ic,jy) - usd(iu,jy))*rdx1     !ix=0,iu=-1,let usd(-1,jy)=usd(1,jy)
        if(dsign(1.1d0,us)>=0.0d0 .and. usd(ix-1,jy)==0.0d0) dudx = 0.0d0
        if(dsign(1.1d0,us)<=0.0d0 .and. usd(ix+1,jy)==0.0d0) dudx = 0.0d0
      
      !.upwinding difference
        vs = 0.5d0*(vnd(ix,jy) + vnd(ix,jy-1))
        ju = jy - int(max(0.,dsign(1.1d0,vs)))
        jc = ju + 1    !jy
        rdy1 = 1.0d0/(ysv(ix,jc)-ysv(ix,ju))
        dudy = (usd(ix,jc) - usd(ix,ju))*rdy1
        if(dsign(1.1d0,vs)>=0.0d0 .and. usd(ix,jy-1)==0.0d0) dudy = 0.0d0
        if(dsign(1.1d0,vs)<=0.0d0 .and. usd(ix,jy+1)==0.0d0) dudy = 0.0d0
      
        ustr(ix,jy) = usd(ix,jy) - dtme*(us*dudx + vs*dudy)

      end do
    end do


    do ix=1,iN
      do jy=0, jN

      !.upwinding difference
        us = 0.5d0*(und(ix,jy)+und(ix-1,jy))
	    iu = ix - int(max(0.,dsign(1.1d0,us)))
        ic = iu + 1
        rdx1 = 1.0d0/(xsh(ic,jy)-xsh(iu,jy))
        phi_x = 1.0d0/(alpha*gamma*xsh(ix,jy)**(gamma-1.0d0))                
        dvdx = phi_x*(vsd(ic,jy) - vsd(iu,jy))*rdx1
        if(dsign(1.1d0,us)>=0.0d0 .and. vsd(ix-1,jy)==0.0d0) dvdx = 0.0d0
        if(dsign(1.1d0,us)<=0.0d0 .and. vsd(ix+1,jy)==0.0d0) dvdx = 0.0d0
        
      !.upwinding difference
        vs = vsd(ix,jy)
        ju = jy - int(max(0.,dsign(1.1d0,vs)))
        jc = ju + 1
        rdy1 = 1.0d0/(ysh(ix,jc)-ysh(ix,ju))
        dvdy = (vsd(ix,jc) - vsd(ix,ju))*rdy1 !jy=0,ju=-1, let vsd(ix,-1)=vsd(ix,1)
        if(dsign(1.1d0,vs)>=0.0d0 .and. vsd(ix,jy-1)==0.0d0) dvdy = 0.0d0
        if(dsign(1.1d0,vs)<=0.0d0 .and. vsd(ix,jy+1)==0.0d0) dvdy = 0.0d0
                
        vstr(ix,jy) = vsd(ix,jy) - dtme*(us*dvdx + vs*dvdy)


      end do
    end do

    return
  
  end subroutine     



  subroutine sbcCDS(dtme)

    implicit real(dbl_kind)(a-h,o-z), integer(i-n)
    
    real(dbl_kind), dimension(:,:), allocatable :: w1,  w2
    
    ustr = 0.0d0;   vstr = 0.0d0
    do jy=1, jN
      do ix=1, iN

        im = ix - 1;   ip = ix + 1
        jm = jy - 1;   jp = jy + 1 
      !.central difference
        r2dx = 1.0d0/(xsv(ip,jy)-xsv(im,jy))
        dudx = (usd(ip,jy) - usd(im,jy))*r2dx
        
        r2dy = 1.0d0/(ysv(ix,jp)-ysv(ix,jm))
        dudy = (usd(ix,jp) - usd(ix,jm))*r2dy

!        us = ucl(ix,jy)
!        vs = vcl(ix,jy)
        us = 0.5d0*(usd(im,jy)+usd(ip,jy))
        vs = 0.25d0*(vsd(ip,jy)+vsd(im,jy)+vsd(ix,jp)+vsd(ix,jm))
        ustr(ix,jy) = usd(ix,jy) - dtme*(us*dudx + vs*dudy)
!if(ustr(ix,jy).GT.1.0E-6) then
!write(6,*) ustr(ix,jy),ix,jy
!end if
      end do
    end do    

    do ix=1, iN
      do jy=1, jN

        im = ix - 1;   ip = ix + 1
        jm = jy - 1;   jp = jy + 1
      !.central difference  
        r2dx = 1.0d0/(xsh(ip,jy)-xsh(im,jy))
        dvdx = (vsd(ip,jy) - vsd(im,jy))*r2dx
       
        r2dy = 1.0d0/(ysh(ix,jp)-ysh(ix,jm))
        dvdy = (vsd(ix,jp) - vsd(ix,jm))*r2dy
        
!        vs = vcl(ix,jy)
!        us = ucl(ix,jy)
        vs = 0.5d0*(vsd(ix,jm)+vsd(ix,jp))
        us = 0.25d0*(usd(ip,jy)+usd(im,jy)+usd(ix,jp)+usd(ix,jm))

        vstr(ix,jy) = vsd(ix,jy) - dtme*(us*dvdx + vs*dvdy)

!        if(vstr(ix,jy) .LE. 1.0E-10)  vstr(ix,jy) = 0.0d0

      end do
    end do
        
    return   
      
  end subroutine sbcCDS
  

  subroutine sbcETA(tme, dtme)

    implicit real(dbl_kind)(a-h,o-z), integer(i-n)
  
  !.solve the continuity for elevations at each cell
  !------------------------------------------------------------------------------------------
  !
  !  jcg jacobi conjugate gradient solver from itpack2d, dsrc2c.f90
  !
  !------------------------------------------------------------------------------------------

    integer,        dimension(12)                 :: iPrm   ! 
    real(dbl_kind),  dimension(12)                 :: rPrm   !  
  
    integer,        dimension(:),   allocatable :: iMtr   ! iMax*jMax + 1
    integer,        dimension(:),   allocatable :: jMtr   ! 5*iMax*jMax
    integer,        dimension(:),   allocatable :: iSpc   ! 3*iMax*jMax  
    real(dbl_kind),  dimension(:),   allocatable :: wSpc   ! 6*nCel+4*MITER
    real(dbl_kind),  dimension(:),   allocatable :: etai   ! nCel
    real(dbl_kind),  dimension(:),   allocatable :: qel2   ! nCel
    real(dbl_kind),  dimension(:),   allocatable :: vMtr   ! 5*nCel

    if(IDB==1) write(21, *) ' call sbcETA...'

    nm = iN*jN   
    nSpc = 6*nm + 4*mItr   !.workspace for JCG

    iStat = 0   
    if(iStat==0) allocate(iMtr(nm+1), stat=iStat)
    if(iStat==0) allocate(jMtr(5*nm), stat=iStat)
    if(iStat==0) allocate(iSpc(3*nm), stat=iStat)
    if(iStat==0) allocate(wSpc(nSpc), stat=iStat)
    if(iStat==0) allocate(etai(nm),   stat=iStat)
    if(iStat==0) allocate(qel2(nm),   stat=iStat)
    if(iStat==0) allocate(vMtr(5*nm), stat=iStat)
    if(iStat/=0) then
      if(IDB==1) write(21, *) ' FATAL: allocation problem for sbcMass!'
      if(IDB==1) write(21, *) ' Program terminates abnormally!'
      stop
    end if

    iMtr = 0;       jMtr = 0;       iSpc = 0
    wSpc = 0.0d0;   vMtr = 0.0d0;   qel2 = 0.0d0;   etai = 0.0d0
  
    call sbcSPM(tme, dtme)

  !.assemble sparse matrix format
    nEQs = 0      !.final index of eqs.
    nNZs = 0      !.number of non-zero entries   
    do jy=1, jN
      do ix=1, iN
        nEQs = nEQs + 1;   nNZs = nNZs + 1   
        iMtr(nEQs) = nNZs    !IA=index of first nonzero in JA
        jMtr(nNZs) = iMap(ix,jy,0)   ! JA=column index of of nonzeros in sparse matrix
        if(iMap(ix,jy,0)/=nEQs) then
          if(IDB==1) write(21, *) ' FATAL: impossible 300', ic
          stop
        end if
        etai(nEQs) = eta(ix,jy)
        qel2(nEQs) = rhs(ix,jy)
        vMtr(nNZs) = spm(ix,jy,0)

        do ij=1, 4
          iEQ = iMap(ix,jy,ij)
          if(iEQ<=nm .and. iEQ>nEQs) then  ! upper triangle
            nNZs = nNZs + 1
            jMtr(nNZs) = iMap(ix,jy,ij)
            vMtr(nNZs) = spm(ix,jy,ij)
          end if
        end do
      end do
    end do    
    iMtr(nEQs+1) = nNZs + 1


  !.call sparse matrix solver
    iPrm = 0;   rPrm = 0.0d0 
    call dfault(iPrm,rPrm)
    iPrm(1)  = mItr   !.maximum iteration
    iPrm(2)  = 0      !.level of output msg
    iPrm(4)  = 33     !.output msg to fort.??
    iPrm(5)  = 0      !.symmetric system
    iPrm(11) = 1      !.no timing
    iPrm(12) = 1      !.error analysis
    rPrm(1)  = zeta   !.stopping criterion
    rPrm(8)  = tolr   !.1.0e-12
 
    iSlv = 1
    select case(iSlv)
      case(1); call jcg(nEQs, iMtr, jMtr, vMtr, qel2, etai, &
                &      iSpc, nSpc, wSpc, iPrm, rPrm, iErr)
      case(2); call jsi(nEQs, iMtr, jMtr, vMtr, qel2, etai, &
                &      iSpc, nSpc, wSpc, iPrm, rPrm, iErr)
      case(3); call ssorcg(nEQs, iMtr, jMtr, vMtr, qel2, etai, &
                &      iSpc, nSpc, wSpc, iPrm, rPrm, iErr)
      case(4); call ssorsi(nEQs, iMtr, jMtr, vMtr, qel2, etai, &
                &      iSpc, nSpc, wSpc, iPrm, rPrm, iErr)
    end select 

  !.update water elevation  
    do jy=1, jN
      do ix=1, iN
        i = iN*(jy-1) + ix
        etan(ix,jy) = etai(i)
!if(etan(ix,jy) .GE. 1.0E-9) then
!write(6,*) etan(ix,jy),ix,jy,0
!end if
      end do
      etan(0,jy) = etan(1,jy); etan(iN+1,jy) = etan(iN,jy)
!      etan(0,jy) = 0.0d0; etan(iN+1,jy) = etan(iN,jy)
    end do

    do ix=0,iN+1
      etan(ix,0) = etan(ix,1)
!      etan(ix,jN+1) = 0.0d0
      etan(ix,jN+1)=etan(ix,jN) ! boundary should at y=0 and y=iN, then this is correct
    end do

    eta = etan
  
    if(allocated(iMtr)) deallocate(iMtr)  
    if(allocated(jMtr)) deallocate(jMtr)  
    if(allocated(iSpc)) deallocate(iSpc)
    if(allocated(wSpc)) deallocate(wSpc)
    if(allocated(vMtr)) deallocate(vMtr) 
    if(allocated(qel2)) deallocate(qel2)
    if(allocated(etai)) deallocate(etai)

    if(IDB==1) write(21, *) ' End Call Mass...'
 
    return
  
  end subroutine sbcETA



  subroutine sbcSPM(tme, dtme)

    use modGLB
    implicit real(dbl_kind)(a-h,o-z), integer(i-n)
  !.generate sparse matrix, and right hand side
    spm = 0.0d0;   rhs = 0.0d0
    do jy=1, jN
      do ix=1, iN
        ip = ix + 1;   im = ix - 1
        jp = jy + 1;   jm = jy - 1
        spm(ix,jy,0) = 1.0d0
        rhs(ix,jy) = eta(ix,jy)
        do ij=1, 4
!          if(iMap(ix,jy,ij)>=1 .and. iMap(ix,jy,ij)<=iN*jN) then
            if(ij==1 .or. ij==3) phi_x0 = 1.0d0/(alpha*gamma*xcl(ix,jy)**(gamma-1.0d0))
            if(ij==2 .or. ij==4) phi_x0 = 1.0d0
            select case(ij)
             case(1) !.ix+1/2, jy
              phi_x = 1.0d0/(alpha*gamma*xsv(ix,jy)**(gamma-1.0d0))       
              rc2c = 1.0d0 / (xcl(ip,jy)-xcl(ix,jy))
              rs2s = 1.0d0 / (xsv(ix,jy)-xsv(im,jy))
              ars = arsv(ix,jy);   gij = xGht(ix,jy)
              dij = dsv(ix,jy);    uij = usd(ix,jy)
              eta_ij = -eta(ip,jy)
             case(2) !.ix, jy+1/2
              phi_x = 1.0d0          
              rc2c = 1.0d0 / (ycl(ix,jp)-ycl(ix,jy))
              rs2s = 1.0d0 / (ysh(ix,jy)-ysh(ix,jm))
              ars = arsh(ix,jy);   gij = yGht(ix,jy)
              dij = dsh(ix,jy);    uij = vsd(ix,jy)
              eta_ij = -eta(ix,jp)
             case(3) !.ix-1/2, jy
              phi_x = 1.0d0 / (alpha*gamma*xsv(im,jy)**(gamma-1.0d0))  
              rc2c = 1.0d0 / (xcl(ix,jy)-xcl(im,jy))    
              rs2s = 1.0d0 / (xsv(ix,jy)-xsv(im,jy))
              ars = arsv(im,jy);   gij = -xGht(im,jy)
              dij = dsv(im,jy);    uij = usd(ix,jy)
              eta_ij = eta(ix,jy)
              case(4) !.ix, jy-1/2
              phi_x = 1.0d0
              rc2c = 1.0d0 / (ycl(ix,jy)-ycl(ix,jm)) 
              rs2s = 1.0d0 / (ysh(ix,jy)-ysh(ix,jm))
              ars = arsh(ix,jm);   gij = -yGht(ix,jm)
              dij = dsh(ix,jm);    uij = vsd(ix,jy)
              eta_ij = eta(ix,jy)
            end select

            if(iMap(ix,jy,ij)>=1 .and. iMap(ix,jy,ij)<=iN*jN) then
            const = phi_x0*phi_x*grv*ars*(dtme*dij)**2.0d0*rc2c*rs2s !dij=h
            spm(ix,jy,0)  = spm(ix,jy,0) + const
            spm(ix,jy,ij) = -const
            end if
            rhs_1 = phi_x0*dtme*rs2s*dij*gij
            rhs(ix,jy) = rhs(ix,jy) - rhs_1 + dtme*rs2s*uij*eta_ij
!          end if
        end do

!        write(6,*) spm(ix,jy,0),ix,jy,0
        if(isnan(spm(ix,jy,0))) then
        write(6,*) rhs(ix,jy),eta_ij,uij,ix,jy, 100
        pause
        end if

      end do
    end do

  !.influence of boundary condition
  !---------------------------------------------------------------------------
  !
  !  1. Periodic BC: cells with periodic BC can be viewed as interior cell, 
  !     which is taken into account in iMap
  !  2. Water Level BC: affects both sparse matrix and R.H.S
  !  3. Velocity BC: only affect R.H.S.
  !  5. Land BC: do not affect the formation of the Matrix
  !  5. Currently, SOUTH and NORTH Boundary can only by Land or Periodic BC
  !
  !---------------------------------------------------------------------------
  !
!    etan = 0.0d0;   usdn = 0.0d0;   vsdn = 0.0d0
!    call sbcBCTS(tme+dtme, dtme)
  !.Boundary condition at WEST, "Dis"
    im = 0;   ix = 1;   ip = 2
    do jy=2, jN
      phi_x = 1.0d0 / (alpha*gamma*xsv(im,jy)**(gamma-1.0d0))
      phi_x0 = 1.0d0/(alpha*gamma*xcl(ix,jy)**(gamma-1.0d0))
      rc2c = 1.0d0 / (xcl(ix,jy)-xcl(im,jy))
      rs2s = 1.0d0 / (xsv(ix,jy)-xsv(im,jy))
      ars = arsv(im,jy);   gij = -xGht(im,jy)
      dij = dsv(im,jy);    uij = usd(ix,jy)
      eta_ij = eta(im,jy);
      select case (cSweBW)
       case('Eta')
       case('Dis')
        const = phi_x0*phi_x*grv*ars*(dtme*dij)**2.0d0*rc2c*rs2s
        spm(ix,jy,0) = spm(ix,jy,0) - const   !coeff of eta(i-1,j) shouldn't included since eta(i-1,j)=eta(i,j)
        spm(ix,jy,3) = 0.0d0
!        rhs_1 = phi_x0*dtme*rs2s*dij*gij
!        rhs(ix,jy) = rhs(ix,jy) - rhs_1 + dtme*rs2s*uij*eta_ij

!        rhs(ix,jy) = rhs(ix,jy) - phi_x0*dtme*dij*usd(ix,jy)*rs2s
       case('Lnd')
       case('Prd')
      end select
    end do
  !.Boundary condition at EAST, "RAD"
    im = iN - 1;   ix = iN;   ip = iN + 1
    do jy=1, jN-1
      phi_x0 = 1.0d0/(alpha*gamma*xcl(ix,jy)**(gamma-1.0d0))
      phi_x = 1.0d0/(alpha*gamma*xsv(ix,jy)**(gamma-1.0d0))
      rc2c = 1.0d0 / (xcl(ip,jy)-xcl(ix,jy))
      rs2s = 1.0d0 / (xsv(ix,jy)-xsv(im,jy))
      ars = arsv(ix,jy);   gij = xGht(ix,jy)
      dij = dsv(ix,jy);    uij = usd(ix,jy)
      eta_ij = eta(ix,jy)
      select case (cSweBE)
       case('Eta')
       case('Dis')
 !       rhs(ix,jy) = rhs(ix,jy) - phi_x0*dtme*dij*usd(ix,jy)*rs2s
       case('Lnd')
       case('Prd')
       case('RAD')
        const = phi_x0*phi_x*grv*ars*(dtme*dij)**2.0d0*rc2c*rs2s
        spm(ix,jy,0) = spm(ix,jy,0) - const
        spm(ix,jy,1) = 0.0d0

!        rhs_1 = phi_x0*dtme*rs2s*dij*gij
!        rhs(ix,jy) = rhs(ix,jy) + dtme*rs2s*uij*eta_ij - rhs_1
      end select
    end do

  !.Boundary condition at South
    jy=1
    do ix=2, iN
      phi_x0 = 1.0d0/(alpha*gamma*xcl(ix,jy)**(gamma-1.0d0))
      phi_x = 1.0d0/(alpha*gamma*xsv(ix,jy)**(gamma-1.0d0))
      rc2c = 1.0d0 / (xcl(ix,jy)-xcl(ix,jm))
      rs2s = 1.0d0 / (xsv(ix,jy)-xsv(ix,jm))
      ars = arsv(ix,jm);   gij = xGht(ix,jy)
      dij = dsv(ix,jm);    uij = usd(ix,jy)
      eta_ij = eta(ix,jy)
      
      select case (cSweBE)
       case('Lnd')  ! already dealed with by iMap

       case('Prd')    !.taken care of by iMap
       end select
    end do
    
 
  !.Boundary condition at North
    jy=jN
    do ix=2, iN-1
      phi_x0 = 1.0d0/(alpha*gamma*xcl(ix,jy)**(gamma-1.0d0))
      phi_x = 1.0d0/(alpha*gamma*xsv(ix,jy)**(gamma-1.0d0))
      rc2c = 1.0d0 / (xcl(ip,jy)-xcl(ix,jy))
      rs2s = 1.0d0 / (xsv(ix,jy)-xsv(ix,jm))
      ars = arsv(ix,jy);   gij = xGht(ix,jy)
      dij = dsv(ix,jy);    uij = usd(ix,jy)
      eta_ij = eta(ix,jy)

      select case (cSweBE)
       case('Lnd')  ! already dealed with by iMap

       case('Prd')    !.taken care of by iMap
       end select
    end do
 
    return

  end subroutine sbcSPM


  subroutine sbcBCUPD(tme, dtme)

    implicit real(dbl_kind)(a-h,o-z), integer(i-n)  
	
    if(IDB==1) write(21,*) ' call sbcBCUPD... '

  !.EAST Boundary
    ix = iN  !move boundary from iN+1 to iN
    do jy=0, jN+1
      select case(cSweBE)
       case('Eta')
       case('Dis')
       case('Lnd')
       case('RAD')
!         eta(ix,jy) = eta(ix-1,jy)
         usd(ix,jy) = usd(ix-1,jy)
         usd(ix+1,jy) = usd(ix,jy) !iN+1 is ghost point
         vsd(ix+1,jy) = vsd(ix,jy)
       case default; write(21,*) " EastBC, No such boundary type:", cSweBE
      end select    
    end do     
  !.WEST Boundary
    ix = 0
    do jy=0, jN+1
      select case(cSweBW)
       case('Eta')
       case('Dis')
        if(ynd(ix,jy) .LE. (ynd(ix,jN/3)-0.026d0)) then
         usd(ix,jy) = 0.013d0
       else  if(ynd(ix,jy) .LE. (ynd(ix,jN/3)+0.012d0)) then
         usd(ix,jy) = 0.013d0+0.037d0*(1.0d0+dtanh((ynd(ix,jy)-ynd(ix,jN/3)+0.01d0)/0.026d0))
        else if(ynd(ix,jy) .LE. (ynd(ix,jN/3)+0.012d0+0.167d0)) then
         l = (ynd(ix,jy)-ynd(ix,jN/3)-0.012d0)/(2.0d0*0.167d0)
         usd(ix,jy) = 0.013d0+0.095d0*(1.5d0*l-0.5d0*l**3)
!write(6,*) l,usd(ix,jy)
        else
         usd(ix,jy) = 0.174d0 
        end if
!        usd(ix,jy) = qxcBW
        usd(-1,jy) = usd(ix,jy) !-1 is ghost point
        vsd(ix,jy) = 0.0d0
        vsd(ix+1,jy) = 0.0d0  !Vx=0 at ix=0
!        eta(ix,jy) = eta(ix+1,jy)
       case('Lnd')
       case default; write(21,*) " WestBC, No such boundary type:", cSweBW
      end select
    end do
  !.NORTH Boundary
    jy = jN  !move the boundary from jN+1 to jN
    do ix=0, iN
      select case(cSweBN)
       case('Eta')
       case('Dis')
       case('Lnd')
!         eta(ix,jy) = eta(ix,jy-1)
         usd(ix,jy+1) = usd(ix,jy)  ! Uy=0 at jy=jN
         vsd(ix,jy) = 0.0d0
         vsd(ix,jy+1) = 0.0d0
       case('Prd')
         eta(ix,jy) = eta(ix,jN)
         usd(ix,jy) = usd(ix,jN)
         vsd(ix,jy) = vsd(ix,jN)
       case default
        write(21,*) " NorthBC, No such boundary type:", cSweBN
      end select
    end do
  !.South Boundary
    jy = 0  
    do ix=0, iN
      select case(cSweBS)
       case('Eta')
       case('Dis')
       case('Lnd')
!         eta(ix,jy) = eta(ix,jy+1)
         usd(ix,jy) = usd(ix,jy+1)
         vsd(ix,jy) = 0.0d0
         vsd(ix,-1) = 0.0d0
       case('Prd')
         eta(ix,jy) = eta(ix,jN)
         usd(ix,jy) = usd(ix,jN)
         vsd(ix,jy) = vsd(ix,jN)
       case default
        write(21,*) " SouthBC, No such boundary type:", cSweBS
      end select
    end do 
    
    if(IDB==1) write(21,*) ' end call sbcBCUPD... '
    
    return
    
  end subroutine sbcBCUPD

end module modSWE
