subroutine SWE() bind (c,name="SWE")
  
  use modGLB;   use modMSH
  use modSWE
  implicit real(dbl_kind)(a-h,o-z), integer(i-n)

!.File Information
  integer,           parameter :: IDCtr = 5
  integer,           parameter :: IDMap = 7
  integer,           parameter :: IDDmp = 8
  integer,           parameter :: IDRST = 9
  integer,           parameter :: IDHot = 10
  character(len=30), parameter :: cfCtr = 'SweCtr.ctr'
  character(len=30), parameter :: cfMap = './result/SweMap.dat'
  character(len=30), parameter :: cfDmp = './result/SweDmp.dmp'
  character(len=30), parameter :: cfRst = './result/SweRst.rst'      
!.Read input data
  integer            :: iHot, iSWE  
  character(len=10)  :: CTRL
  character(len=15)  :: KEYWORD
  character(len=30)  :: cfHot
!.time information
  integer            :: minStp, maxStp 
  real(dbl_kind)      :: tme, dtme, tmeFnl   
  real(dbl_kind)      :: tmeObs, tmeMap
  real(dbl_kind)      :: dtObs, dtMap, dtDmp, dtRst
!. Output data
  real(dbl_kind)      :: stress_secY, P_secY, D_secY, Urms_secY
  real(dbl_kind)      :: stress_042, stress_020
  real(dbl_kind)      :: vort, stress, P, D, Urms
  logical            :: exFile

 !.call system_clock(iStrT, icount_rate)

   inquire(file=cfCtr, exist=exFile)
   if(.not.exFile) then
     write(21,*) 'FATAL: control file does not exist!', cfCtr
     stop
   end if
   
   call constant(pii)

 !.open control file
   iEnd = 0
   open(IDCtr, file=cfCtr, status='OLD')
 !.read from file
   do while(iEnd==0)
   !.reading module parameter  
     read(IDCtr, *, err=910, end=950) CTRL 
	 write(21, *) CTRL 
   !.read for each module 
     select case(CTRL)
      case('MAIN')
     !.include swe module?
       read(IDCtr, *, err=910, end=950) KEYWORD, iSWE
 	   write(21,*) KEYWORD, iSWE
	   if(KEYWORD/='iSWE') goto 930
     !.final simulation time  
       read(IDCtr, *, err=910, end=950) KEYWORD, tmeFnl
	   write(21,*) KEYWORD, tmeFnl
       if(KEYWORD/='tmeFnl') goto 930          
     !.time step
       read(IDCtr, *, err=910, end=950) KEYWORD, dtme
	   write(21,*) KEYWORD, dtme
       if(KEYWORD/='tmeStp') goto 930
     !.observation file time
       read(IDCtr, *, err=910, end=950) KEYWORD, tmeObs
	   write(21,*) KEYWORD, tmeObs
       if(KEYWORD/='tmeObs') goto 930            
       read(IDCtr, *, err=910, end=950) KEYWORD, dtObs
	   write(21,*) KEYWORD, dtObs
       if(KEYWORD/='dtObs') goto 930       
     !.map file time
       read(IDCtr, *, err=910, end=950) KEYWORD, tmeMap
	   write(21,*) KEYWORD, tmeFlc
       if(KEYWORD/='tmeMap') goto 930            
       read(IDCtr, *, err=910, end=950) KEYWORD, dtMap
	   write(21,*) KEYWORD, dtMap
       if(KEYWORD/='dtMap') goto 930            
       read(IDCtr, *, err=910, end=950) KEYWORD, dtDmp
	   write(21,*) KEYWORD, dtDmp
       if(KEYWORD/='dtDmp') goto 930
       read(IDCtr, *, err=910, end=950) KEYWORD, dtRst
	   write(21,*) KEYWORD, dtRst
       if(KEYWORD/='dtRst') goto 930
     !.restart file   
       read(IDCtr, *, err=910, end=950) KEYWORD, iHot
	   write(21,*) KEYWORD, iHot
       if(KEYWORD/='iHot') goto 930
       if(iHot==1) then
         read(IDCtr, *, err=910, end=950) KEYWORD, cfHot
	     write(21,*) KEYWORD, cfHot
         if(KEYWORD/='cfHot') goto 930  
       end if       

      case('GRID') !.pass to grid module
       call sbmRead(IDCtr, cfCtr)
 
      case('SWE') !.pass to shallow water module
       if(iSWE==1) call sbcRead(IDCtr, cfCtr)	  

      case('END') !.input complete, quit file
       write(21,*) ' INFOR: reading complete!'
       iEnd = 1

      case default
       write(21,*) ' Fatal: Control Key not Defined!', CTRL
       stop

     end select

   end do 

   close(IDCtr)

   if(dtMap .ne. 0.0) then 
     open(IDMap, file=cfMap)
     if(IDB==1) open(199, file='WavMap.dat')     
   end if
   if(dtDmp .ne. 0.0) open(IDDmp, file=cfDmp, form='unformatted')
   if(iHot==1) then
     inquire(file=cfHot, exist=exFile)
     if(.not.exFile) then
       write(21,*) ' FATAL: restart file does not exist!', cfHot
       stop
     end if
     open(IDHot, file=cfHot, form='unformatted')
   end if

   if(iSWE/=0) call sbcInit(iHot, IDHot, minStp, tme)

 !  if(iHot==1) then
 !    call sbiSwe2Wav(tme)   
 !    call sbpWav(tme, dtme)
 !    call sbiWav2Swe(tme) 
 !    close(IDHot)
 !  end if
 !  if(iHot==1) then
 !    call sbpEta(tme, dtme)
 !  end if   

 !.output initial obs file 
   if (dtObs/=0.0d0 .and. tme>=tmeObs) then
     if (mod(tme, dtObs) .lt. dtme) then
        call subSctFileIO(130, 'Y', 400, 0)
        call subSctJ(iStp, tme, 400, 130, iSWE)  !ifile=130
!        call subSctI(iStp, tme, 6, 130, iSWE)
     end if
   end if 
 !.output initial map file
   if (dtMap/=0.0d0 .and. tme>=tmeMap) then
     call subMap(minStp, tme, IDMap, iSWE)
   end if
   
   maxStp = int((tmeFnl-tme)/dtme) + minStp
    
   do iStp = minStp, maxStp

     write(6, '(I10, F15.3)') iStp, tme
     write(21,'(I10, F15.3)') iStp, tme

   !.solving shallow water equation 
     if(iSWE/=0 .and. tme>=0.0d0) call sbcSolv(tme, dtme)
     
     tme = tme + dtme      

   !.output obs file 
     if (tme>=tmeObs .and. dtObs/=0.0d0) then
       if (mod(tme, dtObs) .lt. dtme) then
!        call subSctFileIO(130, 'Y', 400, 0)
        call subSctJ(iStp, tme, 400, 130, iSWE)
!        call subSctI(iStp,tme,6,130,iSWE)
       end if
     end if
   !.output map file
     if(tme>=tmeMap) then
     !.output .dat file     
       if (dtMap/=0.0d0 .and. mod(tme, dtMap)<dtme) call subMap(iStp, tme, IDMap, iSWE)
     !.output .dmp file       
     !  if (dtDmp/=0.0d0 .and. mod(tme, dtDmp)<dtme) call subDmp(iStp, tme, IDDmp, iSWE, iWAV)
     !.output .rst file
     !  if (dtRst/=0.0d0 .and. mod(tme, dtRst)<dtme) call subRst(iStp, tme, IDRst, cfRst, iSWE, iWAV)
     end if

   end do   

   close(IDMap)
   close(IDRst)
   
!  call system_clock(iEndT, icount_rate)
!  tmeBtr = real(iEndT - iStrT)/icount_rate
!  write(21, '(A30, f8.2, A12)') '     INFOR: Simulation took', tmeBtr, 'seconds...'

910  write(21,*) ' FATAL: error reading control file,', cfCtr
     stop
930  write(21,*) ' FATAL: KEYWORD does not defined!', KEYWORD
     stop
950  write(21,*) ' FATAL: unexpect end of control file,', cfCtr 
     stop

end subroutine SWE



  subroutine subSctFileIO(ifile, cdir, ij, IOC)
  
    use modGLB;   use modMSH;   use modSWE
    implicit real(dbl_kind)(a-h,o-z), integer(i-n)
    character(len=1) :: cdir
    character(len=40) :: flnm
    character(len=3) :: ftrm   
    
    if(IOC==0) then
      write (ftrm,'(I3.3)') ij
      flnm = './result/'//'sec'//cdir//trim(ftrm)//'_TimeSeries044.dat'
      open(ifile+1, file = flnm)
      flnm = './result/'//'sec'//cdir//trim(ftrm)//'_Stress.dat'
      open(ifile+2, file = flnm)
      flnm = './result/'//'sec'//cdir//trim(ftrm)//'_TimeSeries020.dat'
      open(ifile+3, file = flnm)
    end if
    
    if(IOC==1) then
      close(ifile+1);   close(ifile+2);   close(ifile+3)
    end if

    return
    
  end subroutine



!.write X-sectional file 
  subroutine subSctI(istp, tme, j, ifile, iSWE)

    use modGLB;   use modMSH;   use modSWE
    implicit real(dbl_kind)(a-h,o-z), integer(i-n)
    
      if(iSWE==1) then
!       write(ifile+1,'(F15.1, 600E20.8)') tme, (eta(i,j), i=0, iN+1)
      write(ifile+1,'(F15.1, 600E20.8)') tme, (eta(i,j), i=0, iN+1, 4)
      write(ifile+2,'(F15.1, 600E20.8)') tme, (ucl(i,j), i=0, iN+1, 4)
      write(ifile+3,'(F15.1, 600E20.8)') tme, (vcl(i,j), i=0, iN+1, 4)
      end if

    return
    
  end subroutine subSctI
  
  
!.write Y-sectional file 
  subroutine subSctJ(istp, tme, i, ifile, iSWE)

    use modGLB;   use modMSH;   use modSWE
    implicit real(dbl_kind)(a-h,o-z), integer(i-n)
      
     write(ifile+2,*) 'Variables = "y", "stress_secY", "u", "P_secY", "D_secY"'
      write(ifile+2,*) 'Zone T= "tme- ',tme,'"'

      do j = 1, jN
        do ix = 1, iN
        Um(ix,j) = Um(ix-1,j) + ucl(ix,j); Vm(ix,j) = Vm(ix-1,j) + vcl(ix,j)
        end do
      Um(0,j) = Um(iN,j)/iN; Vm(0,j) = Vm(iN,j)/iN

      stress_secY = -(ucl(i,j)-Um(0,j))*(vcl(i,j)-Vm(0,j))
      rdy2 = 1.0d0/(ycl(i,j+1) - ycl(i,j-1))
      dudy = (usd(i,j+1)-usd(i,j-1))*rdy2
      P_secY = dudy*stress_secY
      Urms_secY = dsqrt(ucl(i,j)**2.0d0+vcl(i,j)**2.0d0)
      if(j .LE. jN/3) then
        D_secY = 0.5d0*25.5d0*Urms_secY**3.0d0
      else
        D_secY = 0.0d0
       end if

      ! this is for sectional at x=800
      write(ifile+2,'(6E15.6,/)') ycl(i,j),stress_secY,ucl(i,j),P_secY,D_secY
      end do

      ! this is for the time series at one point x=800, y=42
      stress_044 = 0.0d0
      stress_020 = 0.0d0
      stress_044 = -(ucl(i,44)-Um(0,44))*(vcl(i,44)-Vm(0,44))
      stress_020 = -(ucl(i,20)-Um(0,20))*(vcl(i,20)-Vm(0,20))
      write(ifile+1,'(F15.1, 6E15.6)') tme, ucl(i,44), vcl(i,44), eta(i,44),stress_044
      write(ifile+3,'(F15.1, 6E15.6)') tme, ucl(i,20), vcl(i,20), eta(i,20),stress_020
    return
    
  end subroutine subSctJ  
  

!.write tecplot mapfile .dat  
  subroutine subMap(istp, tme, ifile, iSWE)
  
    use modGLB;   use modMSH;   use modSWE
    implicit real(dbl_kind)(a-h,o-z), integer(i-n)

    write(ifile,*) 'Variables = "x", "y", "eta", "u", "v","stress","vort","P","D"'
    write(ifile, '(A13, F10.3, A1)') ' Zone T="tme - ', tme, '"'
    write(ifile, '(A3, I5, A6, I5, A11)') ' I=', iN/4, ',   J=', jN/4, ',   f=point'
!    write(ifile,*) 'Zone T="STP:1", STRANDID=1, SOLUTIONTIME=0.01'
!    write(ifile, '(A3, I5, A6, I5, A40)') 'I=', iN, ',J=', jN, ', ZONETYPE=Ordered, DATAPACKING=POINT'
!    write(ifile,*) 'DT=(SINGLE SINGLE)'

    ! average velocity in y direction
    do j=1, jN
      do i=1, iN
        Umean(i,j) = Umean(i-1,j) + ucl(i,j); Vmean(i,j) = Vmean(i-1,j) + vcl(i,j)
      end do
      Umean(0,j) = Umean(iN,j)/iN; Vmean(0,j) = Vmean(iN,j)/iN
    end do

    do j=1, jN, 4
      do i=1, iN, 4
        vort = 0.0d0
        !.vorticity
          rdx2 = 1.0d0/(xcl(i+1,j) - xcl(i-1,j))
          rdy2 = 1.0d0/(ycl(i,j+1) - ycl(i,j-1))
          dudy = (usd(i,j+1)-usd(i,j-1))*rdy2
          dvdx = (vsd(i+1,j)-vsd(i-1,j))*rdx2
          vort = dvdx-dudy
          stress = 0.0d0
          stress = -(ucl(i,j)-Umean(0,j))*(vcl(i,j)-Vmean(0,j))
          P = dudy*stress
          Urms = dsqrt(ucl(i,j)**2.0d0+vcl(i,j)**2.0d0)
          if(j .LE. jN/3) then
          D = 0.5d0*25.5d0*Urms**3.0d0
          else
          D = 0.0d0
          end if 
          write(ifile, '(6E15.6, /)') xcl(i,j), ycl(i,j), eta(i,j), ucl(i,j), vcl(i,j), stress, vort, P, D
          write(ifile,*)
      end do
    end do
  
    return

  end subroutine subMap



!.write binary mapfile .dmp
  subroutine subDmp(istp, tme, ifile, iSWE)
  
    use modGLB;   use modMSH;   use modSWE
    implicit real(dbl_kind)(a-h,o-z), integer(i-n)
    
    if(iSWE/=0) then
       write(ifile) istp, tme, iN, jN	
       write(ifile) ((eta(i,j), i=1, iN), j=1, jN), &
           &         ((usd(i,j), i=1, iN), j=1, jN), &
           &         ((vsd(i,j), i=1, iN), j=1, jN)
    end if
    
    return
    
  end subroutine subDmp
    
    
    
!.write binary restart file .rst
  subroutine subRST(istp, tme, ifile, cfile, iSWE, iWAV)
  
    use modGLB;   use modMSH;   use modSWE
    implicit real(dbl_kind)(a-h,o-z), integer(i-n)
    
    logical           :: exfile
    character(len=30) :: cfile
    
    inquire(file=cfile, exist=exFile)
    if(exFile) then
      open (ifile, file=cfile)    
      close(ifile, status='delete')
    end if
    open(ifile, file=cfile, form='unformatted')    
    if(iSWE/=0) then
       write(ifile) istp, tme, iN, jN	
       write(ifile) ((eta(i,j), i=0, iN+1), j=0, jN+1), &
           &         ((usd(i,j), i=0, iN+1), j=0, jN+1), &
           &         ((vsd(i,j), i=0, iN+1), j=0, jN+1)
    end if
     
    close(ifile)
          
    return
    
   end subroutine subRst
        
