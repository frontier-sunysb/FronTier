      SUBROUTINE JCG (NN,IA,JA,A,RHS,U,IWKSP,NW,WKSP,IPARM,RPARM,IERR)
!       
!     ITPACK 2C MAIN SUBROUTINE  JCG  (JACOBI CONJUGATE GRADIENT)   
!     EACH OF THE MAIN SUBROUTINES:   
!           JCG, JSI, SOR, SSORCG, SSORSI, RSCG, RSSI     
!     CAN BE USED INDEPENDENTLY OF THE OTHERS   
!       
! ... FUNCTION:   
!       
!          THIS SUBROUTINE, JCG, DRIVES THE JACOBI CONJUGATE
!          GRADIENT ALGORITHM.
!       
! ... PARAMETER LIST:       
!       
!          N      INPUT INTEGER.  DIMENSION OF THE MATRIX. (= NN)   
!          IA,JA  INPUT INTEGER VECTORS.  THE TWO INTEGER ARRAYS OF 
!                 THE SPARSE MATRIX REPRESENTATION.       
!          A      INPUT D.P. VECTOR.  THE D.P. ARRAY OF THE SPARSE  
!                 MATRIX REPRESENTATION.
!          RHS    INPUT D.P. VECTOR.  CONTAINS THE RIGHT HAND SIDE  
!                 OF THE MATRIX PROBLEM.
!          U      INPUT/OUTPUT D.P. VECTOR.  ON INPUT, U CONTAINS THE 
!                 INITIAL GUESS TO THE SOLUTION. ON OUTPUT, IT CONTAINS 
!                 THE LATEST ESTIMATE TO THE SOLUTION.    
!          IWKSP  INTEGER VECTOR WORKSPACE OF LENGTH 3*N  
!          NW     INPUT INTEGER.  LENGTH OF AVAILABLE WKSP.  ON OUTPUT, 
!                 IPARM(8) IS AMOUNT USED.      
!          WKSP   D.P. VECTOR USED FOR WORKING SPACE.  JACOBI CONJUGATE 
!                 GRADIENT NEEDS THIS TO BE IN LENGTH AT LEAST      
!                 4*N + 2*ITMAX,  IF ISYM = 0  (SYMMETRIC STORAGE)  
!                 4*N + 4*ITMAX,  IF ISYM = 1  (NONSYMMETRIC STORAGE) 
!                 HERE ITMAX = IPARM(1) AND ISYM = IPARM(5) 
!                 (ITMAX IS THE MAXIMUM ALLOWABLE NUMBER OF ITERATIONS) 
!          IPARM  INTEGER VECTOR OF LENGTH 12.  ALLOWS USER TO SPECIFY
!                 SOME INTEGER PARAMETERS WHICH AFFECT THE METHOD.  
!          RPARM  D.P. VECTOR OF LENGTH 12. ALLOWS USER TO SPECIFY SOME 
!                 D.P. PARAMETERS WHICH AFFECT THE METHOD.
!          IER    OUTPUT INTEGER.  ERROR FLAG. (= IERR)   
!       
! ... JCG SUBPROGRAM REFERENCES:      
!       
!          FROM ITPACK    BISRCH, CHGCON, DETERM, DFAULT, ECHALL,   
!                         ECHOUT, EIGVNS, EIGVSS, EQRT1S, ITERM, TIMER, 
!                         ITJCG, IVFILL, PARCON, PERMAT,  
!                         PERROR5, PERVEC, PJAC, PMULT, PRBNDX,      
!                         PSTOP, QSORT, DAXPY, SBELM, SCAL, DCOPY,  
!                         DDOT, SUM3, UNSCAL, VEVMW, VFILL, VOUT,   
!                         WEVMW, ZBRENT 
!          SYSTEM         DABS, DLOG10, DBLE(AMAX0), DMAX1, MOD, DSQRT
!       
!     VERSION:  ITPACK 2C (MARCH 1982)
!       
!     CODE WRITTEN BY:  DAVID KINCAID, ROGER GRIMES, JOHN RESPESS   
!                       CENTER FOR NUMERICAL ANALYSIS     
!                       UNIVERSITY OF TEXAS     
!                       AUSTIN, TX  78712       
!                       (512) 471-1242
!       
!     FOR ADDITIONAL DETAILS ON THE   
!          (A) SUBROUTINE SEE TOMS ARTICLE 1982 
!          (B) ALGORITHM  SEE CNA REPORT 150    
!       
!     BASED ON THEORY BY:  DAVID YOUNG, DAVID KINCAID, LOU HAGEMAN  
!       
!     REFERENCE THE BOOK:  APPLIED ITERATIVE METHODS      
!                          L. HAGEMAN, D. YOUNG 
!                          ACADEMIC PRESS, 1981 
!       
!     **************************************************  
!     *               IMPORTANT NOTE                   *  
!     *                                                *  
!     *      WHEN INSTALLING ITPACK ROUTINES ON A      *  
!     *  DIFFERENT COMPUTER, RESET SOME OF THE VALUES  *  
!     *  IN  SUBROUTNE DFAULT.   MOST IMPORTANT ARE    *  
!     *                                                *  
!     *   DRELPR      MACHINE RELATIVE PRECISION       *  
!     *   RPARM(1)    STOPPING CRITERION               *  
!     *                                                *  
!     *   ALSO CHANGE SYSTEM-DEPENDENT ROUTINE         *  
!     *   SECOND USED IN TIMER                         *  
!     *                                                *  
!     **************************************************  
!       
!     SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(1),JA(1),IWKSP(1),IPARM(12),NN,NW,IERR   
      DOUBLE PRECISION A(1),RHS(NN),U(NN),WKSP(NW),RPARM(12)
!       
!     SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER IB1,IB2,IB3,IB4,IB5,IDGTS,IER,IERPER,ITMAX1,LOOP,N,NB,N3
      DOUBLE PRECISION DIGIT1,DIGIT2,TEMP,TIME1,TIME2,TOL 
!       
! **** BEGIN: ITPACK COMMON 
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! **** END  : ITPACK COMMON 
!       
! ... VARIABLES IN COMMON BLOCK - ITCOM1
!       
!     IN     - ITERATION NUMBER       
!     IS     - ITERATION NUMBER WHEN PARAMETERS LAST CHANGED
!     ISYM   - SYMMETRIC/NONSYMMETRIC STORAGE FORMAT SWITCH 
!     ITMAX  - MAXIMUM NUMBER OF ITERATIONS ALLOWED       
!     LEVEL  - LEVEL OF OUTPUT CONTROL SWITCH   
!     NOUT   - OUTPUT UNIT NUMBER     
!       
! ... VARIABLES IN COMMON BLOCK - ITCOM2
!       
!     ADAPT  - FULLY ADAPTIVE PROCEDURE SWITCH  
!     BETADT - SWITCH FOR ADAPTIVE DETERMINATION OF BETA  
!     CASEII - ADAPTIVE PROCEDURE CASE SWITCH   
!     HALT   - STOPPING TEST SWITCH   
!     PARTAD - PARTIALLY ADAPTIVE PROCEDURE SWITCH
!       
! ... VARIABLES IN COMMON BLOCK - ITCOM3
!       
!     BDELNM - TWO NORM OF B TIMES DELTA-SUPER-N
!     BETAB  - ESTIMATE FOR THE SPECTRAL RADIUS OF LU MATRIX
!     CME    - ESTIMATE OF LARGEST EIGENVALUE   
!     DELNNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION N      
!     DELSNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION S      
!     FF     - ADAPTIVE PROCEDURE DAMPING FACTOR
!     GAMMA  - ACCELERATION PARAMETER 
!     OMEGA  - OVERRELAXATION PARAMETER FOR SOR AND SSOR  
!     QA     - PSEUDO-RESIDUAL RATIO  
!     QT     - VIRTUAL SPECTRAL RADIUS
!     RHO    - ACCELERATION PARAMETER 
!     RRR    - ADAPTIVE PARAMETER     
!     SIGE   - PARAMETER SIGMA-SUB-E  
!     SME    - ESTIMATE OF SMALLEST EIGENVALUE  
!     SPECR  - SPECTRAL RADIUS ESTIMATE FOR SSOR
!     DRELPR - MACHINE RELATIVE PRECISION       
!     STPTST - STOPPING PARAMETER     
!     UDNM   - TWO NORM OF U
!     ZETA   - STOPPING CRITERION     
!       
! ... INITIALIZE COMMON BLOCKS
!       
      LEVEL = IPARM(2)      
      NOUT = IPARM(4)       
      IF (LEVEL.GE.1) WRITE (NOUT,10) 
   10 FORMAT ('0'///1X,'BEGINNING OF ITPACK SOLUTION MODULE  JCG')  
      IER = 0     
      IF (IPARM(1).LE.0) RETURN       
      N = NN      ! N is matrix
      IF (IPARM(11).EQ.0) TIMJ1 = TIMER(DUMMY)  
      IF (LEVEL.GE.3) GO TO 20
      CALL ECHOUT (IPARM,RPARM,1)     
      GO TO 30    
   20 CALL ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,1) 
   30 TEMP = 5.0D2*DRELPR   
      IF (ZETA.GE.TEMP) GO TO 50      
      IF (LEVEL.GE.1) WRITE (NOUT,40) ZETA,DRELPR,TEMP    
   40 FORMAT ('0','*** W A R N I N G ************'/'0',   &
     &   '    IN ITPACK ROUTINE JCG'/' ','    RPARM(1) =',D10.3,    &
     &   ' (ZETA)'/' ','    A VALUE THIS SMALL MAY HINDER CONVERGENCE '/&
     &   ' ','    SINCE MACHINE PRECISION DRELPR =',D10.3/' ',      &
     &   '    ZETA RESET TO ',D10.3)  
      ZETA = TEMP 
   50 CONTINUE    
      TIME1 = RPARM(9)      
      TIME2 = RPARM(10)     
      DIGIT1 = RPARM(11)    
      DIGIT2 = RPARM(12)    
!       
! ... VERIFY N    
!       
      IF (N.GT.0) GO TO 70  
      IER = 11    
      IF (LEVEL.GE.0) WRITE (NOUT,60) N 
   60 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE JCG '/' ',       &
     &   '    INVALID MATRIX DIMENSION, N =',I8)
      GO TO 370   
   70 CONTINUE    
!       
! ... REMOVE ROWS AND COLUMNS IF REQUESTED      
!       
      IF (IPARM(10).EQ.0) GO TO 90    
      TOL = RPARM(8)
      CALL IVFILL (N,IWKSP,0)   !
      CALL VFILL (N,WKSP,0.0D0)    ! set vector WKSP to 0
      CALL SBELM (N,IA,JA,A,RHS,IWKSP,WKSP,TOL,ISYM,LEVEL,NOUT,IER) 
      IF (IER.EQ.0) GO TO 90
      IF (LEVEL.GE.0) WRITE (NOUT,80) IER,TOL   
   80 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE JCG '/' ',       &
     &   '    ERROR DETECTED IN SUBROUTINE  SBELM '/' ',  &
     &   '    WHICH REMOVES ROWS AND COLUMNS OF SYSTEM '/' ',       &
     &   '    WHEN DIAGONAL ENTRY TOO LARGE  '/' ','    IER = ',I5,5X,&
     &   ' RPARM(8) = ',D10.3,' (TOL)') 
      GO TO 370   
!       
! ... INITIALIZE WKSP BASE ADDRESSES. 
!       
   90 IB1 = 1     
      IB2 = IB1+N 
      IB3 = IB2+N 
      IB4 = IB3+N 
      IB5 = IB4+N 
      IPARM(8) = 4*N+2*ITMAX
      IF (ISYM.NE.0) IPARM(8) = IPARM(8)+2*ITMAX
      IF (NW.GE.IPARM(8)) GO TO 110   
      IER = 12    
      IF (LEVEL.GE.0) WRITE (NOUT,100) NW,IPARM(8)
  100 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE JCG '/' ',       &
     &   '    NOT ENOUGH WORKSPACE AT ',I10/' ','    SET IPARM(8) =',I10&
     &   ,' (NW)')
      GO TO 370   
!       
! ... PERMUTE TO  RED-BLACK SYSTEM IF REQUESTED 
!       
  110 NB = IPARM(9) 
      IF (NB.LT.0) GO TO 170
      N3 = 3*N    
      CALL IVFILL (N3,IWKSP,0)
      CALL PRBNDX (N,NB,IA,JA,IWKSP,IWKSP(IB2),LEVEL,NOUT,IER)      
      IF (IER.EQ.0) GO TO 130 
      IF (LEVEL.GE.0) WRITE (NOUT,120) IER,NB   
  120 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE JCG '/' ',       &
     &   '    ERROR DETECTED IN SUBROUTINE  PRBNDX'/' ',  &
     &   '    WHICH COMPUTES THE RED-BLACK INDEXING'/' ','    IER = ',I5&
     &   ,' IPARM(9) = ',I5,' (NB)')  
      GO TO 370   
!       
! ... PERMUTE MATRIX AND RHS
!       
  130 IF (LEVEL.GE.2) WRITE (NOUT,140) NB       
  140 FORMAT (/10X,'ORDER OF BLACK SUBSYSTEM = ',I5,' (NB)')
      CALL PERMAT (N,IA,JA,A,IWKSP,IWKSP(IB3),ISYM,LEVEL,NOUT,IER)  
      IF (IER.EQ.0) GO TO 160 
      IF (LEVEL.GE.0) WRITE (NOUT,150) IER      
  150 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE JCG '/' ',       &
     &   '    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ',  &
     &   '    WHICH DOES THE RED-BLACK PERMUTATION'/' ','    IER = ',I5)
      GO TO 370   
  160 CALL PERVEC (N,RHS,IWKSP)       
      CALL PERVEC (N,U,IWKSP) 
!       
! ... SCALE LINEAR SYSTEM, U, AND RHS BY THE SQUARE ROOT OF THE     
! ... DIAGONAL ELEMENTS.    
!       
  170 CONTINUE    
      CALL VFILL (IPARM(8),WKSP,0.0D0)
      CALL SCAL (N,IA,JA,A,RHS,U,WKSP,LEVEL,NOUT,IER)  !iER=0, LEVEL=0
      IF (IER.EQ.0) GO TO 190 
      IF (LEVEL.GE.0) WRITE (NOUT,180) IER      
  180 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE JCG '/' ',       &
     &   '    ERROR DETECTED IN SUBROUTINE  SCAL  '/' ',  &
     &   '    WHICH SCALES THE SYSTEM   '/' ','    IER = ',I5)      
      GO TO 370   
  190 IF (LEVEL.LE.2) GO TO 220       
      WRITE (NOUT,200)      
  200 FORMAT (///1X,'IN THE FOLLOWING, RHO AND GAMMA ARE',&
     &   ' ACCELERATION PARAMETERS')  
      IF (ADAPT) WRITE (NOUT,210)     
  210 FORMAT (1X,'CME IS THE ESTIMATE OF THE LARGEST EIGENVALUE OF',&
     &   ' THE JACOBI MATRIX')
  220 IF (IPARM(11).NE.0) GO TO 230   !IPARM(11)=1
      TIMI1 = TIMER(DUMMY)

!       
! ... COMPUTE INITIAL PSEUDO-RESIDUAL 
!       
  230 CONTINUE    
      CALL DCOPY (N,RHS,1,WKSP(IB2),1)
      CALL PJAC (N,IA,JA,A,U,WKSP(IB2))
      CALL VEVMW (N,WKSP(IB2),U)
!       
! ... ITERATION SEQUENCE    
!
      ITMAX1 = ITMAX+1      
      DO 250 LOOP = 1,ITMAX1
         IN = LOOP-1
         IF (MOD(IN,2).EQ.1) GO TO 240
!write(6,*) WKSP(IB1),WKSP(IB2),WKSP(IB3),WKSP(IB4),WKSP(IB5)
!       
! ... CODE FOR THE EVEN ITERATIONS.  IN is odd
!       
!     U           = U(IN)             WKSP(IB2) = DEL(IN) 
!     WKSP(IB1)   = U(IN-1)           WKSP(IB3) = DEL(IN-1) 
!
         CALL ITJCG (N,IA,JA,A,U,WKSP(IB1),WKSP(IB2),WKSP(IB3),WKSP(IB4)&
     &      ,WKSP(IB5))
!       
         IF (HALT) GO TO 280
         GO TO 250
!       
! ... CODE FOR THE ODD ITERATIONS. IN is even
!       
!     U           = U(IN-1)           WKSP(IB2) = DEL(IN-1) 
!     WKSP(IB1)   = U(IN)             WKSP(IB3) = DEL(IN) 
!       
  240    CALL ITJCG (N,IA,JA,A,WKSP(IB1),U,WKSP(IB3),WKSP(IB2),WKSP(IB4)&
     &      ,WKSP(IB5))     
!       
         IF (HALT) GO TO 280

  250 CONTINUE    
!       
! ... ITMAX HAS BEEN REACHED
!       
      IF (IPARM(11).NE.0) GO TO 260   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  260 IER = 13    
      IF (LEVEL.GE.1) WRITE (NOUT,270) ITMAX    
  270 FORMAT ('0','*** W A R N I N G ************'/'0',   &
     &   '    IN ITPACK ROUTINE JCG'/' ','    FAILURE TO CONVERGE IN',I5&
     &   ,' ITERATIONS')    
      IF (IPARM(3).EQ.0) RPARM(1) = STPTST      
      GO TO 310   
!       
! ... METHOD HAS CONVERGED  
!       
  280 IF (IPARM(11).NE.0) GO TO 290   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  290 IF (LEVEL.GE.1) WRITE (NOUT,300) IN       
  300 FORMAT (/1X,'JCG  HAS CONVERGED IN ',I5,' ITERATIONS')
!       
! ... PUT SOLUTION INTO U IF NOT ALREADY THERE. 
!       
  310 CONTINUE    
      IF (MOD(IN,2).EQ.1) CALL DCOPY (N,WKSP(IB1),1,U,1)  
!       
! ... UNSCALE THE MATRIX, SOLUTION, AND RHS VECTORS.      
!       
      CALL UNSCAL (N,IA,JA,A,RHS,U,WKSP)
!       
! ... UN-PERMUTE MATRIX,RHS, AND SOLUTION       
!       
      IF (IPARM(9).LT.0) GO TO 340    
      CALL PERMAT (N,IA,JA,A,IWKSP(IB2),IWKSP(IB3),ISYM,LEVEL,NOUT, &
     &   IERPER)  
      IF (IERPER.EQ.0) GO TO 330      
      IF (LEVEL.GE.0) WRITE (NOUT,320) IERPER   
  320 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE JCG '/' ',       &
     &   '    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ',  &
     &   '    WHICH UNDOES THE RED-BLACK PERMUTATION   '/' ',       &
     &   '    IER = ',I5)   
      IF (IER.EQ.0) IER = IERPER      
      GO TO 370   
  330 CALL PERVEC (N,RHS,IWKSP(IB2))  
      CALL PERVEC (N,U,IWKSP(IB2))    
!       
! ... OPTIONAL ERROR ANALYSIS 
!       
  340 IDGTS = IPARM(12)     
      IF (IDGTS.LT.0) GO TO 350       
      IF (IPARM(2).LE.0) IDGTS = 0    
      CALL PERROR5 (N,IA,JA,A,RHS,U,WKSP,DIGIT1,DIGIT2,IDGTS)
!       
! ... SET RETURN PARAMETERS IN IPARM AND RPARM  
!       
  350 IPARM(8) = IPARM(8)-2*(ITMAX-IN)
      IF (IPARM(11).NE.0) GO TO 360   
      TIMJ2 = TIMER(DUMMY)  
      TIME2 = DBLE(TIMJ2-TIMJ1)       
  360 IF (ISYM.NE.0) IPARM(8) = IPARM(8)-2*(ITMAX-IN)     
      IF (IPARM(3).NE.0) GO TO 370    
      IPARM(1) = IN 
      IPARM(9) = NB 
      RPARM(2) = CME
      RPARM(3) = SME
      RPARM(9) = TIME1      
      RPARM(10) = TIME2     
      RPARM(11) = DIGIT1    
      RPARM(12) = DIGIT2    
!       
  370 CONTINUE    
      IERR = IER  
      IF (LEVEL.GE.3) CALL ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,2)     
!       
      RETURN      
      END








      SUBROUTINE JSI (NN,IA,JA,A,RHS,U,IWKSP,NW,WKSP,IPARM,RPARM,IERR)
!       
!     ITPACK 2C MAIN SUBROUTINE  JSI  (JACOBI SEMI-ITERATIVE)       
!     EACH OF THE MAIN SUBROUTINES:   
!           JCG, JSI, SOR, SSORCG, SSORSI, RSCG, RSSI     
!     CAN BE USED INDEPENDENTLY OF THE OTHERS   
!       
! ... FUNCTION:   
!       
!          THIS SUBROUTINE, JSI, DRIVES THE JACOBI SEMI-  
!          ITERATION ALGORITHM.       
!       
! ... PARAMETER LIST:       
!       
!          N      INPUT INTEGER.  DIMENSION OF THE MATRIX. (= NN)   
!          IA,JA  INPUT INTEGER VECTORS.  THE TWO INTEGER ARRAYS OF 
!                 THE SPARSE MATRIX REPRESENTATION.       
!          A      INPUT D.P. VECTOR.  THE D.P. ARRAY OF THE SPARSE  
!                 MATRIX REPRESENTATION.
!          RHS    INPUT D.P. VECTOR.  CONTAINS THE RIGHT HAND SIDE  
!                 OF THE MATRIX PROBLEM.
!          U      INPUT/OUTPUT D.P. VECTOR.  ON INPUT, U CONTAINS THE 
!                 INITIAL GUESS TO THE SOLUTION. ON OUTPUT, IT CONTAINS 
!                 THE LATEST ESTIMATE TO THE SOLUTION.    
!          IWKSP  INTEGER VECTOR WORKSPACE OF LENGTH 3*N  
!          NW     INPUT INTEGER.  LENGTH OF AVAILABLE WKSP.  ON OUTPUT, 
!                 IPARM(8) IS AMOUNT USED.      
!          WKSP   D.P. VECTOR USED FOR WORKING SPACE.  JACOBI SI    
!                 NEEDS THIS TO BE IN LENGTH AT LEAST     
!                 2*N       
!          IPARM  INTEGER VECTOR OF LENGTH 12.  ALLOWS USER TO SPECIFY
!                 SOME INTEGER PARAMETERS WHICH AFFECT THE METHOD.  
!          RPARM  D.P. VECTOR OF LENGTH 12.  ALLOWS USER TO SPECIFY SOME
!                 D.P. PARAMETERS WHICH AFFECT THE METHOD.
!          IER    OUTPUT INTEGER.  ERROR FLAG. (= IERR)   
!       
! ... JSI SUBPROGRAM REFERENCES:      
!       
!          FROM ITPACK   BISRCH, CHEBY, CHGSI, CHGSME, DFAULT, ECHALL,
!                        ECHOUT, ITERM, TIMER, ITJSI, IVFILL, PAR   
!                        PERMAT, PERROR5, PERVEC, PJAC, PMULT, PRBNDX, 
!                        PSTOP, PVTBV, QSORT, DAXPY, SBELM, SCAL,   
!                        DCOPY, DDOT, SUM3, TSTCHG, UNSCAL, VEVMW,  
!                        VFILL, VOUT, WEVMW     
!          SYSTEM        DABS, DLOG10, DBLE(AMAX0), DMAX1, DBLE(FLOAT), 
!                        MOD,DSQRT    
!       
!     VERSION:  ITPACK 2C (MARCH 1982)
!       
!     CODE WRITTEN BY:  DAVID KINCAID, ROGER GRIMES, JOHN RESPESS   
!                       CENTER FOR NUMERICAL ANALYSIS     
!                       UNIVERSITY OF TEXAS     
!                       AUSTIN, TX  78712       
!                       (512) 471-1242
!       
!     FOR ADDITIONAL DETAILS ON THE   
!          (A) SUBROUTINE SEE TOMS ARTICLE 1982 
!          (B) ALGORITHM  SEE CNA REPORT 150    
!       
!     BASED ON THEORY BY:  DAVID YOUNG, DAVID KINCAID, LOU HAGEMAN  
!       
!     REFERENCE THE BOOK:  APPLIED ITERATIVE METHODS      
!                          L. HAGEMAN, D. YOUNG 
!                          ACADEMIC PRESS, 1981 
!       
!     **************************************************  
!     *               IMPORTANT NOTE                   *  
!     *                                                *  
!     *      WHEN INSTALLING ITPACK ROUTINES ON A      *  
!     *  DIFFERENT COMPUTER, RESET SOME OF THE VALUES  *  
!     *  IN  SUBROUTNE DFAULT.   MOST IMPORTANT ARE    *  
!     *                                                *  
!     *   DRELPR      MACHINE RELATIVE PRECISION       *  
!     *   RPARM(1)    STOPPING CRITERION               *  
!     *                                                *  
!     *   ALSO CHANGE SYSTEM-DEPENDENT ROUTINE         *  
!     *   SECOND USED IN TIMER                         *  
!     *                                                *  
!     **************************************************  
!       
!     SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(1),JA(1),IWKSP(1),IPARM(12),NN,NW,IERR   
      DOUBLE PRECISION A(1),RHS(NN),U(NN),WKSP(NW),RPARM(12)
!       
!     SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER IB1,IB2,IB3,ICNT,IDGTS,IER,IERPER,ITMAX1,LOOP,N,NB,N3 
      DOUBLE PRECISION DIGIT1,DIGIT2,TEMP,TIME1,TIME2,TOL 
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
! ... VARIABLES IN COMMON BLOCK - ITCOM1
!       
!     IN     - ITERATION NUMBER       
!     IS     - ITERATION NUMBER WHEN PARAMETERS LAST CHANGED
!     ISYM   - SYMMETRIC/NONSYMMETRIC STORAGE FORMAT SWITCH 
!     ITMAX  - MAXIMUM NUMBER OF ITERATIONS ALLOWED       
!     LEVEL  - LEVEL OF OUTPUT CONTROL SWITCH   
!     NOUT   - OUTPUT UNIT NUMBER     
!       
! ... VARIABLES IN COMMON BLOCK - ITCOM2
!       
!     ADAPT  - FULLY ADAPTIVE PROCEDURE SWITCH  
!     BETADT - SWITCH FOR ADAPTIVE DETERMINATION OF BETA  
!     CASEII - ADAPTIVE PROCEDURE CASE SWITCH   
!     HALT   - STOPPING TEST SWITCH   
!     PARTAD - PARTIALLY ADAPTIVE PROCEDURE SWITCH
!       
! ... VARIABLES IN COMMON BLOCK - ITCOM3
!       
!     BDELNM - TWO NORM OF B TIMES DELTA-SUPER-N
!     BETAB  - ESTIMATE FOR THE SPECTRAL RADIUS OF LU MATRIX
!     CME    - ESTIMATE OF LARGEST EIGENVALUE   
!     DELNNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION N      
!     DELSNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION S      
!     FF     - ADAPTIVE PROCEDURE DAMPING FACTOR
!     GAMMA  - ACCELERATION PARAMETER 
!     OMEGA  - OVERRELAXATION PARAMETER FOR SOR AND SSOR  
!     QA     - PSEUDO-RESIDUAL RATIO  
!     QT     - VIRTUAL SPECTRAL RADIUS
!     RHO    - ACCELERATION PARAMETER 
!     RRR    - ADAPTIVE PARAMETER     
!     SIGE   - PARAMETER SIGMA-SUB-E  
!     SME    - ESTIMATE OF SMALLEST EIGENVALUE  
!     SPECR  - SPECTRAL RADIUS ESTIMATE FOR SSOR
!     DRELPR - MACHINE RELATIVE PRECISION       
!     STPTST - STOPPING PARAMETER     
!     UDNM   - TWO NORM OF U
!     ZETA   - STOPPING CRITERION     
!       
! ... INITIALIZE COMMON BLOCKS
!       
      LEVEL = IPARM(2)      
      NOUT = IPARM(4)       
      IF (LEVEL.GE.1) WRITE (NOUT,10) 
   10 FORMAT ('0'///1X,'BEGINNING OF ITPACK SOLUTION MODULE  JSI')  
      IER = 0     
      IF (IPARM(1).LE.0) RETURN       
      N = NN      
      IF (IPARM(11).EQ.0) TIMJ1 = TIMER(DUMMY)  
      IF (LEVEL.GE.3) GO TO 20
      CALL ECHOUT (IPARM,RPARM,2)     
      GO TO 30    
   20 CALL ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,1) 
   30 TEMP = 5.0D2*DRELPR   
      IF (ZETA.GE.TEMP) GO TO 50      
      IF (LEVEL.GE.1) WRITE (NOUT,40) ZETA,DRELPR,TEMP    
   40 FORMAT ('0','*** W A R N I N G ************'/'0',   &
     &   '    IN ITPACK ROUTINE JSI'/' ','    RPARM(1) =',D10.3,    &
     &   ' (ZETA)'/' ','    A VALUE THIS SMALL MAY HINDER CONVERGENCE '/&
     &   ' ','    SINCE MACHINE PRECISION DRELPR =',D10.3/' ',      &
     &   '    ZETA RESET TO ',D10.3)  
      ZETA = TEMP 
   50 CONTINUE    
      TIME1 = RPARM(9)      
      TIME2 = RPARM(10)     
      DIGIT1 = RPARM(11)    
      DIGIT2 = RPARM(12)    
!       
! ... VERIFY N    
!       
      IF (N.GT.0) GO TO 70  
      IER = 21    
      IF (LEVEL.GE.0) WRITE (NOUT,60) N 
   60 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE JSI '/' ',       &
     &   '    INVALID MATRIX DIMENSION, N =',I8)
      GO TO 360   
   70 CONTINUE    
!       
! ... REMOVE ROWS AND COLUMNS IF REQUESTED      
!       
      IF (IPARM(10).EQ.0) GO TO 90    
      TOL = RPARM(8)
      CALL IVFILL (N,IWKSP,0) 
      CALL VFILL (N,WKSP,0.0D0)       
      CALL SBELM (N,IA,JA,A,RHS,IWKSP,WKSP,TOL,ISYM,LEVEL,NOUT,IER) 
      IF (IER.EQ.0) GO TO 90
      IF (LEVEL.GE.0) WRITE (NOUT,80) IER,TOL   
   80 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE JSI '/' ',       &
     &   '    ERROR DETECTED IN SUBROUTINE  SBELM '/' ',  &
     &   '    WHICH REMOVES ROWS AND COLUMNS OF SYSTEM '/' ',       &
     &   '    WHEN DIAGONAL ENTRY TOO LARGE  '/' ','    IER = ',I5,5X,&
     &   ' RPARM(8) = ',D10.3,' (TOL)') 
      GO TO 360   
!       
! ... INITIALIZE WKSP BASE ADDRESSES. 
!       
   90 IB1 = 1     
      IB2 = IB1+N 
      IB3 = IB2+N 
      IPARM(8) = 2*N
      IF (NW.GE.IPARM(8)) GO TO 110   
      IER = 22    
      IF (LEVEL.GE.0) WRITE (NOUT,100) NW,IPARM(8)
  100 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE JSI '/' ',       &
     &   '    NOT ENOUGH WORKSPACE AT ',I10/' ','    SET IPARM(8) =',I10&
     &   ,' (NW)')
      GO TO 360   
!       
! ... PERMUTE TO  RED-BLACK SYSTEM IF REQUESTED 
!       
  110 NB = IPARM(9) 
      IF (NB.LT.0) GO TO 170
      N3 = 3*N    
      CALL IVFILL (N3,IWKSP,0)
      CALL PRBNDX (N,NB,IA,JA,IWKSP,IWKSP(IB2),LEVEL,NOUT,IER)      
      IF (IER.EQ.0) GO TO 130 
      IF (LEVEL.GE.0) WRITE (NOUT,120) IER,NB   
  120 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE JSI '/' ',       &
     &   '    ERROR DETECTED IN SUBROUTINE  PRBNDX'/' ',  &
     &   '    WHICH COMPUTES THE RED-BLACK INDEXING'/' ','    IER = ',I5&
     &   ,' IPARM(9) = ',I5,' (NB)')  
      GO TO 360   
!       
! ... PERMUTE MATRIX AND RHS
!       
  130 IF (LEVEL.GE.2) WRITE (NOUT,140) NB       
  140 FORMAT (/10X,'ORDER OF BLACK SUBSYSTEM = ',I5,' (NB)')
      CALL PERMAT (N,IA,JA,A,IWKSP,IWKSP(IB3),ISYM,LEVEL,NOUT,IER)  
      IF (IER.EQ.0) GO TO 160 
      IF (LEVEL.GE.0) WRITE (NOUT,150) IER      
  150 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE JSI '/' ',       &
     &   '    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ',  &
     &   '    WHICH DOES THE RED-BLACK PERMUTATION'/' ','    IER = ',I5)
      GO TO 360   
  160 CALL PERVEC (N,RHS,IWKSP)       
      CALL PERVEC (N,U,IWKSP) 
!       
! ... SCALE LINEAR SYSTEM, U, AND RHS BY THE SQUARE ROOT OF THE     
! ... DIAGONAL ELEMENTS.    
!       
  170 CONTINUE    
      CALL VFILL (IPARM(8),WKSP,0.0D0)
      CALL SCAL (N,IA,JA,A,RHS,U,WKSP,LEVEL,NOUT,IER)     
      IF (IER.EQ.0) GO TO 190 
      IF (LEVEL.GE.0) WRITE (NOUT,180) IER      
  180 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE JSI '/' ',       &
     &   '    ERROR DETECTED IN SUBROUTINE  SCAL  '/' ',  &
     &   '    WHICH SCALES THE SYSTEM   '/' ','    IER = ',I5)      
      GO TO 360   
  190 IF (LEVEL.LE.2) GO TO 210       
      WRITE (NOUT,200)      
  200 FORMAT (///1X,'IN THE FOLLOWING, RHO AND GAMMA ARE',&
     &   ' ACCELERATION PARAMETERS')  
  210 IF (IPARM(11).NE.0) GO TO 220   
      TIMI1 = TIMER(DUMMY)  
!       
! ... ITERATION SEQUENCE    
!       
  220 ITMAX1 = ITMAX+1      
      DO 240 LOOP = 1,ITMAX1
         IN = LOOP-1
         IF (MOD(IN,2).EQ.1) GO TO 230
!       
! ... CODE FOR THE EVEN ITERATIONS.   
!       
!     U           = U(IN)   
!     WKSP(IB1)   = U(IN-1) 
!       
         CALL ITJSI (N,IA,JA,A,RHS,U,WKSP(IB1),WKSP(IB2),ICNT)      
!       
         IF (HALT) GO TO 270
         GO TO 240
!       
! ... CODE FOR THE ODD ITERATIONS.    
!       
!     U           = U(IN-1) 
!     WKSP(IB1)   = U(IN)   
!       
  230    CALL ITJSI (N,IA,JA,A,RHS,WKSP(IB1),U,WKSP(IB2),ICNT)      
!       
         IF (HALT) GO TO 270
  240 CONTINUE    
!       
! ... ITMAX HAS BEEN REACHED
!       
      IF (IPARM(11).NE.0) GO TO 250   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  250 IER = 23    
      IF (LEVEL.GE.1) WRITE (NOUT,260) ITMAX    
  260 FORMAT ('0','*** W A R N I N G ************'/'0',   &
     &   '    IN ITPACK ROUTINE JSI'/' ','    FAILURE TO CONVERGE IN',I5&
     &   ,' ITERATIONS')    
      IF (IPARM(3).EQ.0) RPARM(1) = STPTST      
      GO TO 300   
!       
! ... METHOD HAS CONVERGED  
!       
  270 IF (IPARM(11).NE.0) GO TO 280   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  280 IF (LEVEL.GE.1) WRITE (NOUT,290) IN       
  290 FORMAT (/1X,'JSI  HAS CONVERGED IN ',I5,' ITERATIONS')
!       
! ... PUT SOLUTION INTO U IF NOT ALREADY THERE. 
!       
  300 CONTINUE    
      IF (MOD(IN,2).EQ.1) CALL DCOPY (N,WKSP(IB1),1,U,1)  
!       
! ... UNSCALE THE MATRIX, SOLUTION, AND RHS VECTORS.      
!       
      CALL UNSCAL (N,IA,JA,A,RHS,U,WKSP)
!       
! ... UN-PERMUTE MATRIX,RHS, AND SOLUTION       
!       
      IF (IPARM(9).LT.0) GO TO 330    
      CALL PERMAT (N,IA,JA,A,IWKSP(IB2),IWKSP(IB3),ISYM,LEVEL,NOUT, &
     &   IERPER)  
      IF (IERPER.EQ.0) GO TO 320      
      IF (LEVEL.GE.0) WRITE (NOUT,310) IERPER   
  310 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE JSI '/' ',       &
     &   '    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ',  &
     &   '    WHICH UNDOES THE RED-BLACK PERMUTATION   '/' ',       &
     &   '    IER = ',I5)   
      IF (IER.EQ.0) IER = IERPER      
      GO TO 360   
  320 CALL PERVEC (N,RHS,IWKSP(IB2))  
      CALL PERVEC (N,U,IWKSP(IB2))    
!       
! ... OPTIONAL ERROR ANALYSIS 
!       
  330 IDGTS = IPARM(12)     
      IF (IDGTS.LT.0) GO TO 340       
      IF (IPARM(2).LE.0) IDGTS = 0    
      CALL PERROR5 (N,IA,JA,A,RHS,U,WKSP,DIGIT1,DIGIT2,IDGTS)
!       
! ... SET RETURN PARAMETERS IN IPARM AND RPARM  
!       
  340 IF (IPARM(11).NE.0) GO TO 350   
      TIMJ2 = TIMER(DUMMY)  
      TIME2 = DBLE(TIMJ2-TIMJ1)       
  350 IF (IPARM(3).NE.0) GO TO 360    
      IPARM(1) = IN 
      IPARM(9) = NB 
      RPARM(2) = CME
      RPARM(3) = SME
      RPARM(9) = TIME1      
      RPARM(10) = TIME2     
      RPARM(11) = DIGIT1    
      RPARM(12) = DIGIT2    
!       
  360 CONTINUE    
      IERR = IER  
      IF (LEVEL.GE.3) CALL ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,2)     
!       
      RETURN      
      END 
      SUBROUTINE SOR (NN,IA,JA,A,RHS,U,IWKSP,NW,WKSP,IPARM,RPARM,IERR)
!       
!     ITPACK 2C MAIN SUBROUTINE  SOR  (SUCCESSIVE OVERRELATION)     
!     EACH OF THE MAIN SUBROUTINES:   
!           JCG, JSI, SOR, SSORCG, SSORSI, RSCG, RSSI     
!     CAN BE USED INDEPENDENTLY OF THE OTHERS   
!       
! ... FUNCTION:   
!       
!          THIS SUBROUTINE, SOR, DRIVES THE  SUCCESSIVE   
!          OVERRELAXATION ALGORITHM.  
!       
! ... PARAMETER LIST:       
!       
!          N      INPUT INTEGER.  DIMENSION OF THE MATRIX. (= NN)   
!          IA,JA  INPUT INTEGER VECTORS.  THE TWO INTEGER ARRAYS OF 
!                 THE SPARSE MATRIX REPRESENTATION.       
!          A      INPUT D.P. VECTOR.  THE D.P. ARRAY OF THE SPARSE  
!                 MATRIX REPRESENTATION 
!          RHS    INPUT D.P. VECTOR.  CONTAINS THE RIGHT HAND SIDE  
!                 OF THE MATRIX PROBLEM.
!          U      INPUT/OUTPUT D.P. VECTOR.  ON INPUT, U CONTAINS THE 
!                 INITIAL GUESS TO THE SOLUTION. ON OUTPUT, IT CONTAINS 
!                 THE LATEST ESTIMATE TO THE SOLUTION.    
!          IWKSP  INTEGER VECTOR WORKSPACE OF LENGTH 3*N  
!          NW     INPUT INTEGER.  LENGTH OF AVAILABLE WKSP.  ON OUTPUT, 
!                 IPARM(8) IS AMOUNT USED.      
!          WKSP   D.P. VECTOR USED FOR WORKING SPACE.  SOR NEEDS THIS 
!                 TO BE IN LENGTH AT LEAST  N   
!          IPARM  INTEGER VECTOR OF LENGTH 12.  ALLOWS USER TO SPECIFY
!                 SOME INTEGER PARAMETERS WHICH AFFECT THE METHOD.  
!          RPARM  D.P. VECTOR OF LENGTH 12. ALLOWS USER TO SPECIFY SOME 
!                 D.P. PARAMETERS WHICH AFFECT THE METHOD.
!          IER    OUTPUT INTEGER.  ERROR FLAG. (= IERR)   
!       
! ... SOR SUBPROGRAM REFERENCES:      
!       
!          FROM ITPACK   BISRCH, DFAULT, ECHALL, ECHOUT, IPSTR, ITERM,
!                        TIMER, ITSOR, IVFILL, PERMAT, PERROR5,      
!                        PERVEC, PFSOR1, PMULT, PRBNDX, PSTOP, QSORT, 
!                        SBELM, SCAL, DCOPY, DDOT, TAU, UNSCAL, VFILL,
!                        VOUT, WEVMW  
!          SYSTEM        DABS, DLOG10, DBLE(AMAX0), DMAX1, DBLE(FLOAT), 
!                        DSQRT
!       
!     VERSION:  ITPACK 2C (MARCH 1982)
!       
!     CODE WRITTEN BY:  DAVID KINCAID, ROGER GRIMES, JOHN RESPESS   
!                       CENTER FOR NUMERICAL ANALYSIS     
!                       UNIVERSITY OF TEXAS     
!                       AUSTIN, TX  78712       
!                       (512) 471-1242
!       
!     FOR ADDITIONAL DETAILS ON THE   
!          (A) SUBROUTINE SEE TOMS ARTICLE 1982 
!          (B) ALGORITHM  SEE CNA REPORT 150    
!       
!     BASED ON THEORY BY:  DAVID YOUNG, DAVID KINCAID, LOU HAGEMAN  
!       
!     REFERENCE THE BOOK:  APPLIED ITERATIVE METHODS      
!                          L. HAGEMAN, D. YOUNG 
!                          ACADEMIC PRESS, 1981 
!       
!     **************************************************  
!     *               IMPORTANT NOTE                   *  
!     *                                                *  
!     *      WHEN INSTALLING ITPACK ROUTINES ON A      *  
!     *  DIFFERENT COMPUTER, RESET SOME OF THE VALUES  *  
!     *  IN  SUBROUTNE DFAULT.   MOST IMPORTANT ARE    *  
!     *                                                *  
!     *   DRELPR      MACHINE RELATIVE PRECISION       *  
!     *   RPARM(1)    STOPPING CRITERION               *  
!     *                                                *  
!     *   ALSO CHANGE SYSTEM-DEPENDENT ROUTINE         *  
!     *   SECOND USED IN TIMER                         *  
!     *                                                *  
!     **************************************************  
!       
!     SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(1),JA(1),IWKSP(1),IPARM(12),NN,NW,IERR   
      DOUBLE PRECISION A(1),RHS(NN),U(NN),WKSP(NW),RPARM(12)
!       
!     SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER IB1,IB2,IB3,IDGTS,IER,IERPER,ITMAX1,LOOP,N,NB,N3      
      DOUBLE PRECISION DIGIT1,DIGIT2,TEMP,TIME1,TIME2,TOL 
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
! ... VARIABLES IN COMMON BLOCK - ITCOM1
!       
!     IN     - ITERATION NUMBER       
!     IS     - ITERATION NUMBER WHEN PARAMETERS LAST CHANGED
!     ISYM   - SYMMETRIC/NONSYMMETRIC STORAGE FORMAT SWITCH 
!     ITMAX  - MAXIMUM NUMBER OF ITERATIONS ALLOWED       
!     LEVEL  - LEVEL OF OUTPUT CONTROL SWITCH   
!     NOUT   - OUTPUT UNIT NUMBER     
!       
! ... VARIABLES IN COMMON BLOCK - ITCOM2
!       
!     ADAPT  - FULLY ADAPTIVE PROCEDURE SWITCH  
!     BETADT - SWITCH FOR ADAPTIVE DETERMINATION OF BETA  
!     CASEII - ADAPTIVE PROCEDURE CASE SWITCH   
!     HALT   - STOPPING TEST SWITCH   
!     PARTAD - PARTIALLY ADAPTIVE PROCEDURE SWITCH
!       
! ... VARIABLES IN COMMON BLOCK - ITCOM3
!       
!     BDELNM - TWO NORM OF B TIMES DELTA-SUPER-N
!     BETAB  - ESTIMATE FOR THE SPECTRAL RADIUS OF LU MATRIX
!     CME    - ESTIMATE OF LARGEST EIGENVALUE   
!     DELNNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION N      
!     DELSNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION S      
!     FF     - ADAPTIVE PROCEDURE DAMPING FACTOR
!     GAMMA  - ACCELERATION PARAMETER 
!     OMEGA  - OVERRELAXATION PARAMETER FOR SOR AND SSOR  
!     QA     - PSEUDO-RESIDUAL RATIO  
!     QT     - VIRTUAL SPECTRAL RADIUS
!     RHO    - ACCELERATION PARAMETER 
!     RRR    - ADAPTIVE PARAMETER     
!     SIGE   - PARAMETER SIGMA-SUB-E  
!     SME    - ESTIMATE OF SMALLEST EIGENVALUE  
!     SPECR  - SPECTRAL RADIUS ESTIMATE FOR SSOR
!     DRELPR - MACHINE RELATIVE PRECISION       
!     STPTST - STOPPING PARAMETER     
!     UDNM   - TWO NORM OF U
!     ZETA   - STOPPING CRITERION     
!       
! ... INITIALIZE COMMON BLOCKS
!       
      LEVEL = IPARM(2)      
      NOUT = IPARM(4)       
      IF (LEVEL.GE.1) WRITE (NOUT,10) 
   10 FORMAT ('0'///1X,'BEGINNING OF ITPACK SOLUTION MODULE  SOR')  
      IER = 0     
      IF (IPARM(1).LE.0) RETURN       
      N = NN      
      IF (IPARM(11).EQ.0) TIMJ1 = TIMER(DUMMY)  
      IF (LEVEL.GE.3) GO TO 20
      CALL ECHOUT (IPARM,RPARM,3)     
      GO TO 30    
   20 CALL ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,1) 
   30 TEMP = 5.0D2*DRELPR   
      IF (ZETA.GE.TEMP) GO TO 50      
      IF (LEVEL.GE.1) WRITE (NOUT,40) ZETA,DRELPR,TEMP    
   40 FORMAT ('0','*** W A R N I N G ************'/'0',   &
     &   '    IN ITPACK ROUTINE SOR'/' ','    RPARM(1) =',D10.3,    &
     &   ' (ZETA)'/' ','    A VALUE THIS SMALL MAY HINDER CONVERGENCE '/&
     &   ' ','    SINCE MACHINE PRECISION DRELPR =',D10.3/' ',      &
     &   '    ZETA RESET TO ',D10.3)  
      ZETA = TEMP 
   50 CONTINUE    
      TIME1 = RPARM(9)      
      TIME2 = RPARM(10)     
      DIGIT1 = RPARM(11)    
      DIGIT2 = RPARM(12)    
!       
! ... VERIFY N    
!       
      IF (N.GT.0) GO TO 70  
      IER = 31    
      IF (LEVEL.GE.0) WRITE (NOUT,60) N 
   60 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE SOR '/' ',       &
     &   '    INVALID MATRIX DIMENSION, N =',I8)
      GO TO 360   
   70 CONTINUE    
!       
! ... REMOVE ROWS AND COLUMNS IF REQUESTED      
!       
      IF (IPARM(10).EQ.0) GO TO 90    
      TOL = RPARM(8)
      CALL IVFILL (N,IWKSP,0) 
      CALL VFILL (N,WKSP,0.0D0)       
      CALL SBELM (N,IA,JA,A,RHS,IWKSP,WKSP,TOL,ISYM,LEVEL,NOUT,IER) 
      IF (IER.EQ.0) GO TO 90
      IF (LEVEL.GE.0) WRITE (NOUT,80) IER,TOL   
   80 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE SOR '/' ',       &
     &   '    ERROR DETECTED IN SUBROUTINE  SBELM '/' ',  &
     &   '    WHICH REMOVES ROWS AND COLUMNS OF SYSTEM '/' ',       &
     &   '    WHEN DIAGONAL ENTRY TOO LARGE  '/' ','    IER = ',I5,5X,&
     &   ' RPARM(8) = ',D10.3,' (TOL)') 
      GO TO 360   
!       
! ... INITIALIZE WKSP BASE ADDRESSES. 
!       
   90 IB1 = 1     
      IB2 = IB1+N 
      IB3 = IB2+N 
      IPARM(8) = N
      IF (NW.GE.IPARM(8)) GO TO 110   
      IER = 32    
      IF (LEVEL.GE.0) WRITE (NOUT,100) NW,IPARM(8)
  100 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE SOR '/' ',       &
     &   '    NOT ENOUGH WORKSPACE AT ',I10/' ','    SET IPARM(8) =',I10&
     &   ,' (NW)')
      GO TO 360   
!       
! ... PERMUTE TO  RED-BLACK SYSTEM IF REQUESTED 
!       
  110 NB = IPARM(9) 
      IF (NB.LT.0) GO TO 170
      N3 = 3*N    
      CALL IVFILL (N3,IWKSP,0)
      CALL PRBNDX (N,NB,IA,JA,IWKSP,IWKSP(IB2),LEVEL,NOUT,IER)      
      IF (IER.EQ.0) GO TO 130 
      IF (LEVEL.GE.0) WRITE (NOUT,120) IER,NB   
  120 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE SOR '/' ',       &
     &   '    ERROR DETECTED IN SUBROUTINE  PRBNDX'/' ',  &
     &   '    WHICH COMPUTES THE RED-BLACK INDEXING'/' ','    IER = ',I5&
     &   ,' IPARM(9) = ',I5,' (NB)')  
      GO TO 360   
!       
! ... PERMUTE MATRIX AND RHS
!       
  130 IF (LEVEL.GE.2) WRITE (NOUT,140) NB       
  140 FORMAT (/10X,'ORDER OF BLACK SUBSYSTEM = ',I5,' (NB)')
      CALL PERMAT (N,IA,JA,A,IWKSP,IWKSP(IB3),ISYM,LEVEL,NOUT,IER)  
      IF (IER.EQ.0) GO TO 160 
      IF (LEVEL.GE.0) WRITE (NOUT,150) IER      
  150 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE SOR '/' ',       &
     &   '    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ',  &
     &   '    WHICH DOES THE RED-BLACK PERMUTATION'/' ','    IER = ',I5)
      GO TO 360   
  160 CALL PERVEC (N,RHS,IWKSP)       
      CALL PERVEC (N,U,IWKSP) 
!       
! ... SCALE LINEAR SYSTEM, U, AND RHS BY THE SQUARE ROOT OF THE     
! ... DIAGONAL ELEMENTS.    
!       
  170 CONTINUE    
      CALL VFILL (IPARM(8),WKSP,0.0D0)
      CALL SCAL (N,IA,JA,A,RHS,U,WKSP,LEVEL,NOUT,IER)     
      IF (IER.EQ.0) GO TO 190 
      IF (LEVEL.GE.0) WRITE (NOUT,180) IER      
  180 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE SOR '/' ',       &
     &   '    ERROR DETECTED IN SUBROUTINE  SCAL  '/' ',  &
     &   '    WHICH SCALES THE SYSTEM   '/' ','    IER = ',I5)      
      GO TO 360   
  190 IF (LEVEL.LE.2) GO TO 220       
      IF (ADAPT) WRITE (NOUT,200)     
  200 FORMAT (///1X,'CME IS THE ESTIMATE OF THE LARGEST EIGENVALUE OF', &
     &   ' THE JACOBI MATRIX')
      WRITE (NOUT,210)      
  210 FORMAT (1X,'OMEGA IS THE RELAXATION FACTOR')
  220 IF (IPARM(11).NE.0) GO TO 230   
      TIMI1 = TIMER(DUMMY)  
!       
! ... ITERATION SEQUENCE    
!       
  230 ITMAX1 = ITMAX+1      
      DO 240 LOOP = 1,ITMAX1
         IN = LOOP-1
!       
! ... CODE FOR ONE ITERATION. 
!       
!     U           = U(IN)   
!       
         CALL ITSOR (N,IA,JA,A,RHS,U,WKSP(IB1)) 
!       
         IF (HALT) GO TO 270
  240 CONTINUE    
!       
! ... ITMAX HAS BEEN REACHED
!       
      IF (IPARM(11).NE.0) GO TO 250   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  250 IF (LEVEL.GE.1) WRITE (NOUT,260) ITMAX    
  260 FORMAT ('0','*** W A R N I N G ************'/'0',   &
     &   '    IN ITPACK ROUTINE SOR'/' ','    FAILURE TO CONVERGE IN',I5&
     &   ,' ITERATIONS')    
      IER = 33    
      IF (IPARM(3).EQ.0) RPARM(1) = STPTST      
      GO TO 300   
!       
! ... METHOD HAS CONVERGED  
!       
  270 IF (IPARM(11).NE.0) GO TO 280   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  280 IF (LEVEL.GE.1) WRITE (NOUT,290) IN       
  290 FORMAT (/1X,'SOR  HAS CONVERGED IN ',I5,' ITERATIONS')
!       
! ... UNSCALE THE MATRIX, SOLUTION, AND RHS VECTORS.      
!       
  300 CALL UNSCAL (N,IA,JA,A,RHS,U,WKSP)
!       
! ... UN-PERMUTE MATRIX,RHS, AND SOLUTION       
!       
      IF (IPARM(9).LT.0) GO TO 330    
      CALL PERMAT (N,IA,JA,A,IWKSP(IB2),IWKSP(IB3),ISYM,LEVEL,NOUT, &
     &   IERPER)  
      IF (IERPER.EQ.0) GO TO 320      
      IF (LEVEL.GE.0) WRITE (NOUT,310) IERPER   
  310 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE SOR '/' ',       &
     &   '    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ',  &
     &   '    WHICH UNDOES THE RED-BLACK PERMUTATION   '/' ',       &
     &   '    IER = ',I5)   
      IF (IER.EQ.0) IER = IERPER      
      GO TO 360   
  320 CALL PERVEC (N,RHS,IWKSP(IB2))  
      CALL PERVEC (N,U,IWKSP(IB2))    
!       
! ... OPTIONAL ERROR ANALYSIS 
!       
  330 IDGTS = IPARM(12)     
      IF (IDGTS.LT.0) GO TO 340       
      IF (IPARM(2).LE.0) IDGTS = 0    
      CALL PERROR5 (N,IA,JA,A,RHS,U,WKSP,DIGIT1,DIGIT2,IDGTS)
!       
! ... SET RETURN PARAMETERS IN IPARM AND RPARM  
!       
  340 IF (IPARM(11).NE.0) GO TO 350   
      TIMJ2 = TIMER(DUMMY)  
      TIME2 = DBLE(TIMJ2-TIMJ1)       
  350 IF (IPARM(3).NE.0) GO TO 360    
      IPARM(1) = IN 
      IPARM(9) = NB 
      RPARM(2) = CME
      RPARM(3) = SME
      RPARM(5) = OMEGA      
      RPARM(9) = TIME1      
      RPARM(10) = TIME2     
      RPARM(11) = DIGIT1    
      RPARM(12) = DIGIT2    
!       
  360 CONTINUE    
      IERR = IER  
      IF (LEVEL.GE.3) CALL ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,2)     
!       
      RETURN      
      END 
      SUBROUTINE SSORCG (NN,IA,JA,A,RHS,U,IWKSP,NW,WKSP,IPARM,RPARM,&
     &   IERR)    
!       
!     ITPACK 2C MAIN SUBROUTINE  SSORCG  (SYMMETRIC SUCCESSIVE OVER-
!                                        RELAXATION CONJUGATE GRADIENT) 
!     EACH OF THE MAIN SUBROUTINES:   
!           JCG, JSI, SOR, SSORCG, SSORSI, RSCG, RSSI     
!     CAN BE USED INDEPENDENTLY OF THE OTHERS   
!       
! ... FUNCTION:   
!       
!          THIS SUBROUTINE, SSORCG, DRIVES THE  SYMMETRIC SOR-CG    
!          ALGORITHM.       
!       
! ... PARAMETER LIST:       
!       
!          N      INPUT INTEGER.  DIMENSION OF THE MATRIX. (= NN)   
!          IA,JA  INPUT INTEGER VECTORS.  THE TWO INTEGER ARRAYS OF 
!                 THE SPARSE MATRIX REPRESENTATION.       
!          A      INPUT D.P. VECTOR.  THE D.P. ARRAY OF THE SPARSE  
!                 MATRIX REPRESENTATION.
!          RHS    INPUT D.P. VECTOR.  CONTAINS THE RIGHT HAND SIDE  
!                 OF THE MATRIX PROBLEM.
!          U      INPUT/OUTPUT D.P. VECTOR.  ON INPUT, U CONTAINS THE 
!                 INITIAL GUESS TO THE SOLUTION. ON OUTPUT, IT CONTAINS 
!                 THE LATEST ESTIMATE TO THE SOLUTION.    
!          IWKSP  INTEGER VECTOR WORKSPACE OF LENGTH 3*N  
!          NW     INPUT INTEGER.  LENGTH OF AVAILABLE WKSP.  ON OUTPUT, 
!                 IPARM(8) IS AMOUNT USED.      
!          WKSP   D.P. VECTOR USED FOR WORKING SPACE.  SSOR-CG      
!                 NEEDS TO BE IN LENGTH AT LEAST
!                 6*N + 2*ITMAX,  IF IPARM(5)=0  (SYMMETRIC STORAGE)
!                 6*N + 4*ITMAX,  IF IPARM(5)=1  (NONSYMMETRIC STORAGE) 
!          IPARM  INTEGER VECTOR OF LENGTH 12.  ALLOWS USER TO SPECIFY
!                 SOME INTEGER PARAMETERS WHICH AFFECT THE METHOD.  IF
!          RPARM  D.P. VECTOR OF LENGTH 12. ALLOWS USER TO SPECIFY SOME 
!                 D.P. PARAMETERS WHICH AFFECT THE METHOD.
!          IER    OUTPUT INTEGER.  ERROR FLAG. (= IERR)   
!       
! ... SSORCG SUBPROGRAM REFERENCES:   
!       
!          FROM ITPACK    BISRCH, CHGCON, DETERM, DFAULT, ECHALL,   
!                         ECHOUT, EIGVNS, EIGVSS, EQRT1S, ITERM, TIMER, 
!                         ITSRCG, IVFILL, OMEG, OMGCHG, OMGSTR,     
!                         PARCON, PBETA, PBSOR, PERMAT, PERROR5,     
!                         PERVEC, PFSOR, PJAC, PMULT, PRBNDX, PSTOP, PVT
!                         QSORT, SBELM, SCAL, DCOPY, DDOT, SUM3,    
!                         UNSCAL, VEVMW, VEVPW, VFILL, VOUT, WEVMW, 
!                         ZBRENT      
!          SYSTEM         DABS, DLOG, DLOG10, DBLE(AMAX0), DMAX1, AMIN1,
!                         MOD, DSQRT  
!       
!     VERSION:  ITPACK 2C (MARCH 1982)
!       
!     CODE WRITTEN BY:  DAVID KINCAID, ROGER GRIMES, JOHN RESPESS   
!                       CENTER FOR NUMERICAL ANALYSIS     
!                       UNIVERSITY OF TEXAS     
!                       AUSTIN, TX  78712       
!                       (512) 471-1242
!       
!     FOR ADDITIONAL DETAILS ON THE   
!          (A) SUBROUTINE SEE TOMS ARTICLE 1982 
!          (B) ALGORITHM  SEE CNA REPORT 150    
!       
!     BASED ON THEORY BY:  DAVID YOUNG, DAVID KINCAID, LOU HAGEMAN  
!       
!     REFERENCE THE BOOK:  APPLIED ITERATIVE METHODS      
!                          L. HAGEMAN, D. YOUNG 
!                          ACADEMIC PRESS, 1981 
!       
!     **************************************************  
!     *               IMPORTANT NOTE                   *  
!     *                                                *  
!     *      WHEN INSTALLING ITPACK ROUTINES ON A      *  
!     *  DIFFERENT COMPUTER, RESET SOME OF THE VALUES  *  
!     *  IN  SUBROUTNE DFAULT.   MOST IMPORTANT ARE    *  
!     *                                                *  
!     *   DRELPR      MACHINE RELATIVE PRECISION       *  
!     *   RPARM(1)    STOPPING CRITERION               *  
!     *                                                *  
!     *   ALSO CHANGE SYSTEM-DEPENDENT ROUTINE         *  
!     *   SECOND USED IN TIMER                         *  
!     *                                                *  
!     **************************************************  
!       
!     SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(1),JA(1),IWKSP(1),IPARM(12),NN,NW,IERR   
      DOUBLE PRECISION A(1),RHS(NN),U(NN),WKSP(NW),RPARM(12)
!       
!     SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER IB1,IB2,IB3,IB4,IB5,IB6,IB7,IDGTS,IER,IERPER,ITMAX1,LOOP,N&
     &   ,NB,N3   
      DOUBLE PRECISION BETNEW,DIGIT1,DIGIT2,PBETA,TEMP,TIME1,TIME2,TOL
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
! ... VARIABLES IN COMMON BLOCK - ITCOM1
!       
!     IN     - ITERATION NUMBER       
!     IS     - ITERATION NUMBER WHEN PARAMETERS LAST CHANGED
!     ISYM   - SYMMETRIC/NONSYMMETRIC STORAGE FORMAT SWITCH 
!     ITMAX  - MAXIMUM NUMBER OF ITERATIONS ALLOWED       
!     LEVEL  - LEVEL OF OUTPUT CONTROL SWITCH   
!     NOUT   - OUTPUT UNIT NUMBER     
!       
! ... VARIABLES IN COMMON BLOCK - ITCOM2
!       
!     ADAPT  - FULLY ADAPTIVE PROCEDURE SWITCH  
!     BETADT - SWITCH FOR ADAPTIVE DETERMINATION OF BETA  
!     CASEII - ADAPTIVE PROCEDURE CASE SWITCH   
!     HALT   - STOPPING TEST SWITCH   
!     PARTAD - PARTIALLY ADAPTIVE PROCEDURE SWITCH
!       
! ... VARIABLES IN COMMON BLOCK - ITCOM3
!       
!     BDELNM - TWO NORM OF B TIMES DELTA-SUPER-N
!     BETAB  - ESTIMATE FOR THE SPECTRAL RADIUS OF LU MATRIX
!     CME    - ESTIMATE OF LARGEST EIGENVALUE   
!     DELNNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION N      
!     DELSNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION S      
!     FF     - ADAPTIVE PROCEDURE DAMPING FACTOR
!     GAMMA  - ACCELERATION PARAMETER 
!     OMEGA  - OVERRELAXATION PARAMETER FOR SOR AND SSOR  
!     QA     - PSEUDO-RESIDUAL RATIO  
!     QT     - VIRTUAL SPECTRAL RADIUS
!     RHO    - ACCELERATION PARAMETER 
!     RRR    - ADAPTIVE PARAMETER     
!     SIGE   - PARAMETER SIGMA-SUB-E  
!     SME    - ESTIMATE OF SMALLEST EIGENVALUE  
!     SPECR  - SPECTRAL RADIUS ESTIMATE FOR SSOR
!     DRELPR - MACHINE RELATIVE PRECISION       
!     STPTST - STOPPING PARAMETER     
!     UDNM   - TWO NORM OF U
!     ZETA   - STOPPING CRITERION     
!       
! ... INITIALIZE COMMON BLOCKS
!       
      LEVEL = IPARM(2)      
      NOUT = IPARM(4)       
      IF (IPARM(9).GE.0) IPARM(6) = 2 
      IF (LEVEL.GE.1) WRITE (NOUT,10) 
   10 FORMAT ('0'///1X,'BEGINNING OF ITPACK SOLUTION MODULE  SSORCG') 
      IER = 0     
      IF (IPARM(1).LE.0) RETURN       
      N = NN      
      IF (IPARM(11).EQ.0) TIMJ1 = TIMER(DUMMY)  
      IF (LEVEL.GE.3) GO TO 20
      CALL ECHOUT (IPARM,RPARM,4)     
      GO TO 30    
   20 CALL ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,1) 
   30 TEMP = 5.0D2*DRELPR   
      IF (ZETA.GE.TEMP) GO TO 50      
      IF (LEVEL.GE.1) WRITE (NOUT,40) ZETA,DRELPR,TEMP    
   40 FORMAT ('0','*** W A R N I N G ************'/'0',   &
     &   '    IN ITPACK ROUTINE SSORCG'/' ','    RPARM(1) =',D10.3, &
     &   ' (ZETA)'/' ','    A VALUE THIS SMALL MAY HINDER CONVERGENCE '/&
     &   ' ','    SINCE MACHINE PRECISION DRELPR =',D10.3/' ',      &
     &   '    ZETA RESET TO ',D10.3)  
      ZETA = TEMP 
   50 CONTINUE    
      TIME1 = RPARM(9)      
      TIME2 = RPARM(10)     
      DIGIT1 = RPARM(11)    
      DIGIT2 = RPARM(12)    
!       
! ... VERIFY N    
!       
      IF (N.GT.0) GO TO 70  
      IER = 41    
      IF (LEVEL.GE.0) WRITE (NOUT,60) N 
   60 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE SSORCG '/' ',    &
     &   '    INVALID MATRIX DIMENSION, N =',I8)
      GO TO 390   
   70 CONTINUE    
!       
! ... REMOVE ROWS AND COLUMNS IF REQUESTED      
!       
      IF (IPARM(10).EQ.0) GO TO 90    
      TOL = RPARM(8)
      CALL IVFILL (N,IWKSP,0) 
      CALL VFILL (N,WKSP,0.0D0)       
      CALL SBELM (N,IA,JA,A,RHS,IWKSP,WKSP,TOL,ISYM,LEVEL,NOUT,IER) 
      IF (IER.EQ.0) GO TO 90
      IF (LEVEL.GE.0) WRITE (NOUT,80) IER,TOL   
   80 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE SSORCG '/' ',    &
     &   '    ERROR DETECTED IN SUBROUTINE  SBELM '/' ',  &
     &   '    WHICH REMOVES ROWS AND COLUMNS OF SYSTEM '/' ',       &
     &   '    WHEN DIAGONAL ENTRY TOO LARGE  '/' ','    IER = ',I5,5X,&
     &   ' RPARM(8) = ',D10.3,' (TOL)') 
      GO TO 390   
!       
! ... INITIALIZE WKSP BASE ADDRESSES. 
!       
   90 IB1 = 1     
      IB2 = IB1+N 
      IB3 = IB2+N 
      IB4 = IB3+N 
      IB5 = IB4+N 
      IB6 = IB5+N 
      IB7 = IB6+N 
      IPARM(8) = 6*N+2*ITMAX
      IF (ISYM.NE.0) IPARM(8) = IPARM(8)+2*ITMAX
      IF (NW.GE.IPARM(8)) GO TO 110   
      IER = 42    
      IF (LEVEL.GE.0) WRITE (NOUT,100) NW,IPARM(8)
  100 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE SSORCG '/' ',    &
     &   '    NOT ENOUGH WORKSPACE AT ',I10/' ','    SET IPARM(8) =',I10&
     &   ,' (NW)')
      GO TO 390   
!       
! ... PERMUTE TO  RED-BLACK SYSTEM IF REQUESTED 
!       
  110 NB = IPARM(9) 
      IF (NB.LT.0) GO TO 170
      N3 = 3*N    
      CALL IVFILL (N3,IWKSP,0)
      CALL PRBNDX (N,NB,IA,JA,IWKSP,IWKSP(IB2),LEVEL,NOUT,IER)      
      IF (IER.EQ.0) GO TO 130 
      IF (LEVEL.GE.0) WRITE (NOUT,120) IER,NB   
  120 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE SSORCG '/' ',    &
     &   '    ERROR DETECTED IN SUBROUTINE  PRBNDX'/' ',  &
     &   '    WHICH COMPUTES THE RED-BLACK INDEXING'/' ','    IER = ',I5&
     &   ,' IPARM(9) = ',I5,' (NB)')  
      GO TO 390   
!       
! ... PERMUTE MATRIX AND RHS
!       
  130 IF (LEVEL.GE.2) WRITE (NOUT,140) NB       
  140 FORMAT (/10X,'ORDER OF BLACK SUBSYSTEM = ',I5,' (NB)')
      CALL PERMAT (N,IA,JA,A,IWKSP,IWKSP(IB3),ISYM,LEVEL,NOUT,IER)  
      IF (IER.EQ.0) GO TO 160 
      IF (LEVEL.GE.0) WRITE (NOUT,150) IER      
  150 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE SSORCG '/' ',    &
     &   '    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ',  &
     &   '    WHICH DOES THE RED-BLACK PERMUTATION'/' ','    IER = ',I5)
      GO TO 390   
  160 CALL PERVEC (N,RHS,IWKSP)       
      CALL PERVEC (N,U,IWKSP) 
!       
! ... SCALE LINEAR SYSTEM, U, AND RHS BY THE SQUARE ROOT OF THE     
! ... DIAGONAL ELEMENTS.    
!       
  170 CONTINUE    
      CALL VFILL (IPARM(8),WKSP,0.0D0)
      CALL SCAL (N,IA,JA,A,RHS,U,WKSP,LEVEL,NOUT,IER)     
      IF (IER.EQ.0) GO TO 190 
      IF (LEVEL.GE.0) WRITE (NOUT,180) IER      
  180 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE SSORCG '/' ',    &
     &   '    ERROR DETECTED IN SUBROUTINE  SCAL  '/' ',  &
     &   '    WHICH SCALES THE SYSTEM   '/' ','    IER = ',I5)      
      GO TO 390   
  190 IF (LEVEL.LE.2) GO TO 220       
      WRITE (NOUT,200)      
  200 FORMAT (///1X,'IN THE FOLLOWING, RHO AND GAMMA ARE',&
     &   ' ACCELERATION PARAMETERS')  
      WRITE (NOUT,210)      
  210 FORMAT (1X,'S-PRIME IS AN INITIAL ESTIMATE FOR NEW CME')      
  220 IF (IPARM(11).NE.0) GO TO 230   
      TIMI1 = TIMER(DUMMY)  
!       
! ... SPECIAL PROCEDURE FOR FULLY ADAPTIVE CASE.
!       
  230 CONTINUE    
      IF (.NOT.ADAPT) GO TO 250       
      IF (.NOT.BETADT) GO TO 240      
      CALL VFILL (N,WKSP(IB1),1.D0)   
      BETNEW = PBETA(N,IA,JA,A,WKSP(IB1),WKSP(IB2),WKSP(IB3))/      &
     &   DBLE(FLOAT(N))     
      BETAB = DMAX1(BETAB,.25D0,BETNEW) 
  240 CALL OMEG (0.D0,1)    
      IS = 0      
!       
! ... INITIALIZE FORWARD PSEUDO-RESIDUAL
!       
  250 CALL DCOPY (N,RHS,1,WKSP(IB1),1)
      CALL DCOPY (N,U,1,WKSP(IB2),1)  
      CALL PFSOR (N,IA,JA,A,WKSP(IB2),WKSP(IB1))
      CALL VEVMW (N,WKSP(IB2),U)      
!       
! ... ITERATION SEQUENCE    
!       
      ITMAX1 = ITMAX+1      
      DO 270 LOOP = 1,ITMAX1
         IN = LOOP-1
         IF (MOD(IN,2).EQ.1) GO TO 260
!       
! ... CODE FOR THE EVEN ITERATIONS.   
!       
!     U           = U(IN)       WKSP(IB2) = C(IN) 
!     WKSP(IB1)   = U(IN-1)     WKSP(IB3) = C(IN-1)       
!       
         CALL ITSRCG (N,IA,JA,A,RHS,U,WKSP(IB1),WKSP(IB2),WKSP(IB3),&
     &      WKSP(IB4),WKSP(IB5),WKSP(IB6),WKSP(IB7))      
!       
         IF (HALT) GO TO 300
         GO TO 270
!       
! ... CODE FOR THE ODD ITERATIONS.    
!       
!     U           = U(IN-1)     WKSP(IB2) = C(IN-1)       
!     WKSP(IB1)   = U(IN)       WKSP(IB3) =C(IN)
!       
  260    CALL ITSRCG (N,IA,JA,A,RHS,WKSP(IB1),U,WKSP(IB3),WKSP(IB2),&
     &      WKSP(IB4),WKSP(IB5),WKSP(IB6),WKSP(IB7))      
!       
         IF (HALT) GO TO 300
  270 CONTINUE    
!       
! ... ITMAX HAS BEEN REACHED
!       
      IF (IPARM(11).NE.0) GO TO 280   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  280 IF (LEVEL.GE.1) WRITE (NOUT,290) ITMAX    
  290 FORMAT ('0','*** W A R N I N G ************'/'0',   &
     &   '    IN ITPACK ROUTINE SSORCG'/' ','    FAILURE TO CONVERGE IN'&
     &   ,I5,' ITERATIONS') 
      IER = 43    
      IF (IPARM(3).EQ.0) RPARM(1) = STPTST      
      GO TO 330   
!       
! ... METHOD HAS CONVERGED  
!       
  300 IF (IPARM(11).NE.0) GO TO 310   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  310 IF (LEVEL.GE.1) WRITE (NOUT,320) IN       
  320 FORMAT (/1X,'SSORCG  HAS CONVERGED IN ',I5,' ITERATIONS')     
!       
! ... PUT SOLUTION INTO U IF NOT ALREADY THERE. 
!       
  330 CONTINUE    
      IF (MOD(IN,2).EQ.1) CALL DCOPY (N,WKSP(IB1),1,U,1)  
!       
! ... UNSCALE THE MATRIX, SOLUTION, AND RHS VECTORS.      
!       
      CALL UNSCAL (N,IA,JA,A,RHS,U,WKSP)
!       
! ... UN-PERMUTE MATRIX,RHS, AND SOLUTION       
!       
      IF (IPARM(9).LT.0) GO TO 360    
      CALL PERMAT (N,IA,JA,A,IWKSP(IB2),IWKSP(IB3),ISYM,LEVEL,NOUT, &
     &   IERPER)  
      IF (IERPER.EQ.0) GO TO 350      
      IF (LEVEL.GE.0) WRITE (NOUT,340) IERPER   
  340 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE SSORCG '/' ',    &
     &   '    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ',  &
     &   '    WHICH UNDOES THE RED-BLACK PERMUTATION   '/' ',       &
     &   '    IER = ',I5)   
      IF (IER.EQ.0) IER = IERPER      
      GO TO 390   
  350 CALL PERVEC (N,RHS,IWKSP(IB2))  
      CALL PERVEC (N,U,IWKSP(IB2))    
!       
! ... OPTIONAL ERROR ANALYSIS 
!       
  360 IDGTS = IPARM(12)     
      IF (IDGTS.LT.0) GO TO 370       
      IF (IPARM(2).LE.0) IDGTS = 0    
      CALL PERROR5 (N,IA,JA,A,RHS,U,WKSP,DIGIT1,DIGIT2,IDGTS)
!       
! ... SET RETURN PARAMETERS IN IPARM AND RPARM  
!       
  370 IF (IPARM(11).NE.0) GO TO 380   
      TIMJ2 = TIMER(DUMMY)  
      TIME2 = DBLE(TIMJ2-TIMJ1)       
  380 IPARM(8) = IPARM(8)-2*(ITMAX-IN)
      IF (ISYM.NE.0) IPARM(8) = IPARM(8)-2*(ITMAX-IN)     
      IF (IPARM(3).NE.0) GO TO 390    
      IPARM(1) = IN 
      IPARM(9) = NB 
      RPARM(2) = CME
      RPARM(3) = SME
      RPARM(5) = OMEGA      
      RPARM(6) = SPECR      
      RPARM(7) = BETAB      
      RPARM(9) = TIME1      
      RPARM(10) = TIME2     
      RPARM(11) = DIGIT1    
      RPARM(12) = DIGIT2    
!       
  390 CONTINUE    
      IERR = IER  
      IF (LEVEL.GE.3) CALL ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,2)     
!       
      RETURN      
      END 
      SUBROUTINE SSORSI (NN,IA,JA,A,RHS,U,IWKSP,NW,WKSP,IPARM,RPARM,&
     &   IERR)    
!       
!     ITPACK 2C MAIN SUBROUTINE  SSORSI  (SYMMETRIC SUCCESSIVE RELAX- 
!                                         ATION SEMI-ITERATION)     
!     EACH OF THE MAIN SUBROUTINES:   
!           JCG, JSI, SOR, SSORCG, SSORSI, RSCG, RSSI     
!     CAN BE USED INDEPENDENTLY OF THE OTHERS   
!       
! ... FUNCTION:   
!       
!          THIS SUBROUTINE, SSORSI, DRIVES THE  SYMMETRIC SOR-SI    
!          ALGORITHM.       
!       
! ... PARAMETER LIST:       
!       
!          N      INPUT INTEGER.  DIMENSION OF THE MATRIX. (= NN)   
!          IA,JA  INPUT INTEGER VECTORS.  THE TWO INTEGER ARRAYS OF 
!                 THE SPARSE MATRIX REPRESENTATION.       
!          A      INPUT D.P. VECTOR.  THE D.P. ARRAY OF THE SPARSE  
!                 MATRIX REPRESENTATION.
!          RHS    INPUT D.P. VECTOR.  CONTAINS THE RIGHT HAND SIDE  
!                 OF THE MATRIX PROBLEM.
!          U      INPUT/OUTPUT D.P. VECTOR.  ON INPUT, U CONTAINS THE 
!                 INITIAL GUESS TO THE SOLUTION. ON OUTPUT, IT CONTAINS 
!                 THE LATEST ESTIMATE TO THE SOLUTION.    
!          IWKSP  INTEGER VECTOR WORKSPACE OF LENGTH 3*N  
!          NW     INPUT INTEGER.  LENGTH OF AVAILABLE WKSP.  ON OUTPUT, 
!                 IPARM(8) IS AMOUNT USED.      
!          WKSP   D.P. VECTOR USED FOR WORKING SPACE.  SSORSI       
!                 NEEDS THIS TO BE IN LENGTH AT LEAST  5*N
!          IPARM  INTEGER VECTOR OF LENGTH 12.  ALLOWS USER TO SPECIFY
!                 SOME INTEGER PARAMETERS WHICH AFFECT THE METHOD.  IF
!          RPARM  D.P. VECTOR OF LENGTH 12. ALLOWS USER TO SPECIFY SOME 
!                 D.P. PARAMETERS WHICH AFFECT THE METHOD.
!          IER    OUTPUT INTEGER.  ERROR FLAG. (= IERR)   
!       
! ... SSORSI SUBPROGRAM REFERENCES:   
!       
!          FROM ITPACK    BISRCH, CHEBY, CHGSI, DFAULT, ECHALL, ECHOUT, 
!                         ITERM, TIMER, ITSRSI, IVFILL, OMEG,       
!                         OMGSTR, PARSI, PBETA, PERMAT, PERROR5,     
!                         PERVEC, PFSOR, PMULT, PRBNDX, PSSOR1,     
!                         PSTOP, PVTBV, QSORT, SBELM, SCAL, DCOPY,  
!                         DDOT, SUM3, TSTCHG, UNSCAL, VEVPW, VFILL, 
!                         VOUT, WEVMW 
!          SYSTEM         DABS, DLOG, DLOG10, DBLE(AMAX0), DMAX1, DBLE(F
!                         MOD, DSQRT  
!       
!     VERSION:  ITPACK 2C (MARCH 1982)
!       
!     CODE WRITTEN BY:  DAVID KINCAID, ROGER GRIMES, JOHN RESPESS   
!                       CENTER FOR NUMERICAL ANALYSIS     
!                       UNIVERSITY OF TEXAS     
!                       AUSTIN, TX  78712       
!                       (512) 471-1242
!       
!     FOR ADDITIONAL DETAILS ON THE   
!          (A) SUBROUTINE SEE TOMS ARTICLE 1982 
!          (B) ALGORITHM  SEE CNA REPORT 150    
!       
!     BASED ON THEORY BY:  DAVID YOUNG, DAVID KINCAID, LOU HAGEMAN  
!       
!     REFERENCE THE BOOK:  APPLIED ITERATIVE METHODS      
!                          L. HAGEMAN, D. YOUNG 
!                          ACADEMIC PRESS, 1981 
!       
!     **************************************************  
!     *               IMPORTANT NOTE                   *  
!     *                                                *  
!     *      WHEN INSTALLING ITPACK ROUTINES ON A      *  
!     *  DIFFERENT COMPUTER, RESET SOME OF THE VALUES  *  
!     *  IN  SUBROUTNE DFAULT.   MOST IMPORTANT ARE    *  
!     *                                                *  
!     *   DRELPR      MACHINE RELATIVE PRECISION       *  
!     *   RPARM(1)    STOPPING CRITERION               *  
!     *                                                *  
!     *   ALSO CHANGE SYSTEM-DEPENDENT ROUTINE         *  
!     *   SECOND USED IN TIMER                         *  
!     *                                                *  
!     **************************************************  
!       
!     SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(1),JA(1),IWKSP(1),IPARM(12),NN,NW,IERR   
      DOUBLE PRECISION A(1),RHS(NN),U(NN),WKSP(NW),RPARM(12)
!       
!     SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER IB1,IB2,IB3,IB4,IB5,IDGTS,IER,IERPER,ITMAX1,LOOP,N,NB,N3
      DOUBLE PRECISION BETNEW,DIGIT1,DIGIT2,PBETA,TEMP,TIME1,TIME2,TOL
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
! ... VARIABLES IN COMMON BLOCK - ITCOM1
!       
!     IN     - ITERATION NUMBER       
!     ISYM   - SYMMETRIC/NONSYMMETRIC STORAGE FORMAT SWITCH 
!     IS     - ITERATION NUMBER WHEN PARAMETERS LAST CHANGED
!     ITMAX  - MAXIMUM NUMBER OF ITERATIONS ALLOWED       
!     LEVEL  - LEVEL OF OUTPUT CONTROL SWITCH   
!     NOUT   - OUTPUT UNIT NUMBER     
!       
! ... VARIABLES IN COMMON BLOCK - ITCOM2
!       
!     ADAPT  - FULLY ADAPTIVE PROCEDURE SWITCH  
!     BETADT - SWITCH FOR ADAPTIVE DETERMINATION OF BETA  
!     CASEII - ADAPTIVE PROCEDURE CASE SWITCH   
!     HALT   - STOPPING TEST SWITCH   
!     PARTAD - PARTIALLY ADAPTIVE PROCEDURE SWITCH
!       
! ... VARIABLES IN COMMON BLOCK - ITCOM3
!       
!     BDELNM - TWO NORM OF B TIMES DELTA-SUPER-N
!     BETAB  - ESTIMATE FOR THE SPECTRAL RADIUS OF LU MATRIX
!     CME    - ESTIMATE OF LARGEST EIGENVALUE   
!     DELNNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION N      
!     DELSNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION S      
!     FF     - ADAPTIVE PROCEDURE DAMPING FACTOR
!     GAMMA  - ACCELERATION PARAMETER 
!     OMEGA  - OVERRELAXATION PARAMETER FOR SOR AND SSOR  
!     QA     - PSEUDO-RESIDUAL RATIO  
!     QT     - VIRTUAL SPECTRAL RADIUS
!     RHO    - ACCELERATION PARAMETER 
!     RRR    - ADAPTIVE PARAMETER     
!     SIGE   - PARAMETER SIGMA-SUB-E  
!     SME    - ESTIMATE OF SMALLEST EIGENVALUE  
!     SPECR  - SPECTRAL RADIUS ESTIMATE FOR SSOR
!     DRELPR - MACHINE RELATIVE PRECISION       
!     STPTST - STOPPING PARAMETER     
!     UDNM   - TWO NORM OF U
!     ZETA   - STOPPING CRITERION     
!       
! ... INITIALIZE COMMON BLOCKS
!       
      LEVEL = IPARM(2)      
      NOUT = IPARM(4)       
      IF (IPARM(9).GE.0) IPARM(6) = 2 
      IF (LEVEL.GE.1) WRITE (NOUT,10) 
   10 FORMAT ('0'///1X,'BEGINNING OF ITPACK SOLUTION MODULE  SSORSI') 
      IER = 0     
      IF (IPARM(1).LE.0) RETURN       
      N = NN      
      IF (IPARM(11).EQ.0) TIMJ1 = TIMER(DUMMY)  
      IF (LEVEL.GE.3) GO TO 20
      CALL ECHOUT (IPARM,RPARM,5)     
      GO TO 30    
   20 CALL ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,1) 
   30 TEMP = 5.0D2*DRELPR   
      IF (ZETA.GE.TEMP) GO TO 50      
      IF (LEVEL.GE.1) WRITE (NOUT,40) ZETA,DRELPR,TEMP    
   40 FORMAT ('0','*** W A R N I N G ************'/'0',   &
     &   '    IN ITPACK ROUTINE SSORSI'/' ','    RPARM(1) =',D10.3, &
     &   ' (ZETA)'/' ','    A VALUE THIS SMALL MAY HINDER CONVERGENCE '/&
     &   ' ','    SINCE MACHINE PRECISION DRELPR =',D10.3/' ',      &
     &   '    ZETA RESET TO ',D10.3)  
      ZETA = TEMP 
   50 CONTINUE    
      TIME1 = RPARM(9)      
      TIME2 = RPARM(10)     
      DIGIT1 = RPARM(11)    
      DIGIT2 = RPARM(12)    
!       
! ... VERIFY N    
!       
      IF (N.GT.0) GO TO 70  
      IER = 51    
      IF (LEVEL.GE.0) WRITE (NOUT,60) N 
   60 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE SSORSI '/' ',    &
     &   '    INVALID MATRIX DIMENSION, N =',I8)
      GO TO 380   
   70 CONTINUE    
!       
! ... REMOVE ROWS AND COLUMNS IF REQUESTED      
!       
      IF (IPARM(10).EQ.0) GO TO 90    
      TOL = RPARM(8)
      CALL IVFILL (N,IWKSP,0) 
      CALL VFILL (N,WKSP,0.0D0)       
      CALL SBELM (N,IA,JA,A,RHS,IWKSP,WKSP,TOL,ISYM,LEVEL,NOUT,IER) 
      IF (IER.EQ.0) GO TO 90
      IF (LEVEL.GE.0) WRITE (NOUT,80) IER,TOL   
   80 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE SSORSI '/' ',    &
     &   '    ERROR DETECTED IN SUBROUTINE  SBELM '/' ',  &
     &   '    WHICH REMOVES ROWS AND COLUMNS OF SYSTEM '/' ',       &
     &   '    WHEN DIAGONAL ENTRY TOO LARGE  '/' ','    IER = ',I5,5X,&
     &   ' RPARM(8) = ',D10.3,' (TOL)') 
      GO TO 380   
!       
! ... INITIALIZE WKSP BASE ADDRESSES. 
!       
   90 IB1 = 1     
      IB2 = IB1+N 
      IB3 = IB2+N 
      IB4 = IB3+N 
      IB5 = IB4+N 
      IPARM(8) = 5*N
      IF (NW.GE.IPARM(8)) GO TO 110   
      IER = 52    
      IF (LEVEL.GE.0) WRITE (NOUT,100) NW,IPARM(8)
  100 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE SSORSI '/' ',    &
     &   '    NOT ENOUGH WORKSPACE AT ',I10/' ','    SET IPARM(8) =',I10&
     &   ,' (NW)')
!       
! ... PERMUTE TO  RED-BLACK SYSTEM IF REQUESTED 
!       
  110 NB = IPARM(9) 
      IF (NB.LT.0) GO TO 170
      N3 = 3*N    
      CALL IVFILL (N3,IWKSP,0)
      CALL PRBNDX (N,NB,IA,JA,IWKSP,IWKSP(IB2),LEVEL,NOUT,IER)      
      IF (IER.EQ.0) GO TO 130 
      IF (LEVEL.GE.0) WRITE (NOUT,120) IER,NB   
  120 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE SSORSI '/' ',    &
     &   '    ERROR DETECTED IN SUBROUTINE  PRBNDX'/' ',  &
     &   '    WHICH COMPUTES THE RED-BLACK INDEXING'/' ','    IER = ',I5&
     &   ,' IPARM(9) = ',I5,' (NB)')  
      GO TO 380   
!       
! ... PERMUTE MATRIX AND RHS
!       
  130 IF (LEVEL.GE.2) WRITE (NOUT,140) NB       
  140 FORMAT (/10X,'ORDER OF BLACK SUBSYSTEM = ',I5,' (NB)')
      CALL PERMAT (N,IA,JA,A,IWKSP,IWKSP(IB3),ISYM,LEVEL,NOUT,IER)  
      IF (IER.EQ.0) GO TO 160 
      IF (LEVEL.GE.0) WRITE (NOUT,150) IER      
  150 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE SSORSI '/' ',    &
     &   '    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ',  &
     &   '    WHICH DOES THE RED-BLACK PERMUTATION'/' ','    IER = ',I5)
      GO TO 380   
  160 CALL PERVEC (N,RHS,IWKSP)       
      CALL PERVEC (N,U,IWKSP) 
!       
! ... SCALE LINEAR SYSTEM, U, AND RHS BY THE SQUARE ROOT OF THE     
! ... DIAGONAL ELEMENTS.    
!       
  170 CONTINUE    
      CALL VFILL (IPARM(8),WKSP,0.0D0)
      CALL SCAL (N,IA,JA,A,RHS,U,WKSP,LEVEL,NOUT,IER)     
      IF (IER.EQ.0) GO TO 190 
      IF (LEVEL.GE.0) WRITE (NOUT,180) IER      
  180 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE SSORSI '/' ',    &
     &   '    ERROR DETECTED IN SUBROUTINE  SCAL  '/' ',  &
     &   '    WHICH SCALES THE SYSTEM   '/' ','    IER = ',I5)      
      GO TO 380   
  190 IF (LEVEL.LE.2) GO TO 210       
      WRITE (NOUT,200)      
  200 FORMAT (///1X,'IN THE FOLLOWING, RHO AND GAMMA ARE',&
     &   ' ACCELERATION PARAMETERS')  
  210 IF (IPARM(11).NE.0) GO TO 220   
      TIMI1 = TIMER(DUMMY)  
!       
! ... SPECIAL PROCEDURE FOR FULLY ADAPTIVE CASE.
!       
  220 CONTINUE    
      IF (.NOT.ADAPT) GO TO 240       
      IF (.NOT.BETADT) GO TO 230      
      CALL VFILL (N,WKSP(IB1),1.D0)   
      BETNEW = PBETA(N,IA,JA,A,WKSP(IB1),WKSP(IB2),WKSP(IB3))/      &
     &   DBLE(FLOAT(N))     
      BETAB = DMAX1(BETAB,.25D0,BETNEW) 
  230 CALL OMEG (0.D0,1)    
      IS = 0      
!       
! ... ITERATION SEQUENCE    
!       
  240 ITMAX1 = ITMAX+1      
      DO 260 LOOP = 1,ITMAX1
         IN = LOOP-1
         IF (MOD(IN,2).EQ.1) GO TO 250
!       
! ... CODE FOR THE EVEN ITERATIONS.   
!       
!     U           = U(IN)   
!     WKSP(IB1)   = U(IN-1) 
!       
         CALL ITSRSI (N,IA,JA,A,RHS,U,WKSP(IB1),WKSP(IB2),WKSP(IB3),&
     &      WKSP(IB4),WKSP(IB5))      
!       
         IF (HALT) GO TO 290
         GO TO 260
!       
! ... CODE FOR THE ODD ITERATIONS.    
!       
!     U           = U(IN-1) 
!     WKSP(IB1)   = U(IN)   
!       
  250    CALL ITSRSI (N,IA,JA,A,RHS,WKSP(IB1),U,WKSP(IB2),WKSP(IB3),&
     &      WKSP(IB4),WKSP(IB5))      
!       
         IF (HALT) GO TO 290
  260 CONTINUE    
!       
! ... ITMAX HAS BEEN REACHED
!       
      IF (IPARM(11).NE.0) GO TO 270   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  270 IF (LEVEL.GE.1) WRITE (NOUT,280) ITMAX    
  280 FORMAT ('0','*** W A R N I N G ************'/'0',   &
     &   '    IN ITPACK ROUTINE SSORSI'/' ','    FAILURE TO CONVERGE IN'&
     &   ,I5,' ITERATIONS') 
      IER = 53    
      IF (IPARM(3).EQ.0) RPARM(1) = STPTST      
      GO TO 320   
!       
! ... METHOD HAS CONVERGED  
!       
  290 IF (IPARM(11).NE.0) GO TO 300   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  300 IF (LEVEL.GE.1) WRITE (NOUT,310) IN       
  310 FORMAT (/1X,'SSORSI  HAS CONVERGED IN ',I5,' ITERATIONS')     
!       
! ... PUT SOLUTION INTO U IF NOT ALREADY THERE. 
!       
  320 CONTINUE    
      IF (MOD(IN,2).EQ.1) CALL DCOPY (N,WKSP(IB1),1,U,1)  
!       
! ... UNSCALE THE MATRIX, SOLUTION, AND RHS VECTORS.      
!       
      CALL UNSCAL (N,IA,JA,A,RHS,U,WKSP)
!       
! ... UN-PERMUTE MATRIX,RHS, AND SOLUTION       
!       
      IF (IPARM(9).LT.0) GO TO 350    
      CALL PERMAT (N,IA,JA,A,IWKSP(IB2),IWKSP(IB3),ISYM,LEVEL,NOUT, &
     &   IERPER)  
      IF (IERPER.EQ.0) GO TO 340      
      IF (LEVEL.GE.0) WRITE (NOUT,330) IERPER   
  330 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE SSORSI '/' ',    &
     &   '    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ',  &
     &   '    WHICH UNDOES THE RED-BLACK PERMUTATION   '/' ',       &
     &   '    IER = ',I5)   
      IF (IER.EQ.0) IER = IERPER      
      GO TO 380   
  340 CALL PERVEC (N,RHS,IWKSP(IB2))  
      CALL PERVEC (N,U,IWKSP(IB2))    
!       
! ... OPTIONAL ERROR ANALYSIS 
!       
  350 IDGTS = IPARM(12)     
      IF (IDGTS.LT.0) GO TO 360       
      IF (IPARM(2).LE.0) IDGTS = 0    
      CALL PERROR5 (N,IA,JA,A,RHS,U,WKSP,DIGIT1,DIGIT2,IDGTS)
!       
! ... SET RETURN PARAMETERS IN IPARM AND RPARM  
!       
  360 IF (IPARM(11).NE.0) GO TO 370   
      TIMJ2 = TIMER(DUMMY)  
      TIME2 = DBLE(TIMJ2-TIMJ1)       
  370 IF (IPARM(3).NE.0) GO TO 380    
      IPARM(1) = IN 
      IPARM(9) = NB 
      RPARM(2) = CME
      RPARM(3) = SME
      RPARM(5) = OMEGA      
      RPARM(6) = SPECR      
      RPARM(7) = BETAB      
      RPARM(9) = TIME1      
      RPARM(10) = TIME2     
      RPARM(11) = DIGIT1    
      RPARM(12) = DIGIT2    
!       
  380 CONTINUE    
      IERR = IER  
      IF (LEVEL.GE.3) CALL ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,2)     
!       
      RETURN      
      END 
      SUBROUTINE RSCG (NN,IA,JA,A,RHS,U,IWKSP,NW,WKSP,IPARM,RPARM,IERR) 
!       
!     ITPACK 2C MAIN SUBROUTINE  RSCG  (REDUCED SYSTEM CONJUGATE    
!                                       GRADIENT) 
!     EACH OF THE MAIN SUBROUTINES:   
!           JCG, JSI, SOR, SSORCG, SSORSI, RSCG, RSSI     
!     CAN BE USED INDEPENDENTLY OF THE OTHERS   
!       
! ... FUNCTION:   
!       
!          THIS SUBROUTINE, RSCG, DRIVES THE  REDUCED SYSTEM CG     
!          ALGORITHM.       
!       
! ... PARAMETER LIST:       
!       
!          N     INPUT INTEGER.  DIMENSION OF THE MATRIX. (= NN)    
!                 IN THE RED-BLACK MATRIX.      
!          IA,JA  INPUT INTEGER VECTORS.  THE TWO INTEGER ARRAYS OF 
!                 THE SPARSE MATRIX REPRESENTATION.       
!          A      INPUT D.P. VECTOR.  THE D.P. ARRAY OF THE SPARSE  
!                 MATRIX REPRESENTATION.
!          RHS    INPUT D.P. VECTOR.  CONTAINS THE RIGHT HAND SIDE  
!                 OF THE MATRIX PROBLEM.
!          U      INPUT/OUTPUT D.P. VECTOR.  ON INPUT, U CONTAINS THE 
!                 INITIAL GUESS TO THE SOLUTION. ON OUTPUT, IT CONTAINS 
!                 THE LATEST ESTIMATE TO THE SOLUTION.    
!          IWKSP  INTEGER VECTOR WORKSPACE OF LENGTH 3*N  
!          NW     INPUT INTEGER.  LENGTH OF AVAILABLE WKSP.  ON OUTPUT, 
!                 IPARM(8) IS AMOUNT USED.      
!          WKSP   D.P. VECTOR USED FOR WORKING SPACE.  RSCG NEEDS   
!                 THIS TO BE IN LENGTH AT LEAST 
!                 N+3*NB+2*ITMAX, IF IPARM(5)=0  (SYMMETRIC STORAGE)
!                 N+3*NB+4*ITMAX, IF IPARM(5)=1  (NONSYMMETRIC STORAGE) 
!                 HERE NB IS THE ORDER OF THE BLACK SUBSYSTEM       
!          IPARM  INTEGER VECTOR OF LENGTH 12.  ALLOWS USER TO SPECIFY
!                 SOME INTEGER PARAMETERS WHICH AFFECT THE METHOD.  IF
!          RPARM  D.P. VECTOR OF LENGTH 12. ALLOWS USER TO SPECIFY SOME 
!                 D.P. PARAMETERS WHICH AFFECT THE METHOD.
!          IER    OUTPUT INTEGER. ERROR FLAG. (= IERR)    
!       
! ... RSCG SUBPROGRAM REFERENCES:     
!       
!          FROM ITPACK    BISRCH, CHGCON, DETERM, DFAULT, ECHALL,   
!                         ECHOUT, EIGVNS, EIGVSS, EQRT1S, ITERM, TIMER
!                         ITRSCG, IVFILL, PARCON, PERMAT, 
!                         PERROR5, PERVEC, PMULT, PRBNDX, PRSBLK,    
!                         PRSRED, PSTOP, QSORT, SBELM, SCAL, DCOPY, 
!                         DDOT, SUM3, UNSCAL, VFILL, VOUT, WEVMW,   
!                         ZBRENT      
!          SYSTEM         DABS, DLOG10, DBLE(AMAX0), DMAX1, MOD, DSQRT
!       
!     VERSION:  ITPACK 2C (MARCH 1982)
!       
!     CODE WRITTEN BY:  DAVID KINCAID, ROGER GRIMES, JOHN RESPESS   
!                       CENTER FOR NUMERICAL ANALYSIS     
!                       UNIVERSITY OF TEXAS     
!                       AUSTIN, TX  78712       
!                       (512) 471-1242
!       
!     FOR ADDITIONAL DETAILS ON THE   
!          (A) SUBROUTINE SEE TOMS ARTICLE 1982 
!          (B) ALGORITHM  SEE CNA REPORT 150    
!       
!     BASED ON THEORY BY:  DAVID YOUNG, DAVID KINCAID, LOU HAGEMAN  
!       
!     REFERENCE THE BOOK:  APPLIED ITERATIVE METHODS      
!                          L. HAGEMAN, D. YOUNG 
!                          ACADEMIC PRESS, 1981 
!       
!     **************************************************  
!     *               IMPORTANT NOTE                   *  
!     *                                                *  
!     *      WHEN INSTALLING ITPACK ROUTINES ON A      *  
!     *  DIFFERENT COMPUTER, RESET SOME OF THE VALUES  *  
!     *  IN  SUBROUTNE DFAULT.   MOST IMPORTANT ARE    *  
!     *                                                *  
!     *   DRELPR      MACHINE RELATIVE PRECISION       *  
!     *   RPARM(1)    STOPPING CRITERION               *  
!     *                                                *  
!     *   ALSO CHANGE SYSTEM-DEPENDENT ROUTINE         *  
!     *   SECOND USED IN TIMER                         *  
!     *                                                *  
!     **************************************************  
!       
!     SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(1),JA(1),IWKSP(1),IPARM(12),NN,NW,IERR   
      DOUBLE PRECISION A(1),RHS(NN),U(NN),WKSP(NW),RPARM(12)
!       
!     SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER IB1,IB2,IB3,IB4,IB5,IDGTS,IER,IERPER,ITMAX1,JB3,LOOP,N,NB,&
     &   NR,NRP1,N3 
      DOUBLE PRECISION DIGIT1,DIGIT2,TEMP,TIME1,TIME2,TOL 
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
! ... VARIABLES IN COMMON BLOCK - ITCOM1
!       
!     IN     - ITERATION NUMBER       
!     IS     - ITERATION NUMBER WHEN PARAMETERS LAST CHANGED
!     ISYM   - SYMMETRIC/NONSYMMETRIC STORAGE FORMAT SWITCH 
!     ITMAX  - MAXIMUM NUMBER OF ITERATIONS ALLOWED       
!     LEVEL  - LEVEL OF OUTPUT CONTROL SWITCH   
!     NOUT   - OUTPUT UNIT NUMBER     
!       
! ... VARIABLES IN COMMON BLOCK - ITCOM2
!       
!     ADAPT  - FULLY ADAPTIVE PROCEDURE SWITCH  
!     BETADT - SWITCH FOR ADAPTIVE DETERMINATION OF BETA  
!     CASEII - ADAPTIVE PROCEDURE CASE SWITCH   
!     HALT   - STOPPING TEST SWITCH   
!     PARTAD - PARTIALLY ADAPTIVE PROCEDURE SWITCH
!       
! ... VARIABLES IN COMMON BLOCK - ITCOM3
!       
!     BDELNM - TWO NORM OF B TIMES DELTA-SUPER-N
!     BETAB  - ESTIMATE FOR THE SPECTRAL RADIUS OF LU MATRIX
!     CME    - ESTIMATE OF LARGEST EIGENVALUE   
!     DELNNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION N      
!     DELSNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION S      
!     FF     - ADAPTIVE PROCEDURE DAMPING FACTOR
!     GAMMA  - ACCELERATION PARAMETER 
!     OMEGA  - OVERRELAXATION PARAMETER FOR SOR AND SSOR  
!     QA     - PSEUDO-RESIDUAL RATIO  
!     QT     - VIRTUAL SPECTRAL RADIUS
!     RHO    - ACCELERATION PARAMETER 
!     RRR    - ADAPTIVE PARAMETER     
!     SIGE   - PARAMETER SIGMA-SUB-E  
!     SME    - ESTIMATE OF SMALLEST EIGENVALUE  
!     SPECR  - SPECTRAL RADIUS ESTIMATE FOR SSOR
!     DRELPR - MACHINE RELATIVE PRECISION       
!     STPTST - STOPPING PARAMETER     
!     UDNM   - TWO NORM OF U
!     ZETA   - STOPPING CRITERION     
!       
! ... INITIALIZE COMMON BLOCKS
!       
      LEVEL = IPARM(2)      
      NOUT = IPARM(4)       
      IF (LEVEL.GE.1) WRITE (NOUT,10) 
   10 FORMAT ('0'///1X,'BEGINNING OF ITPACK SOLUTION MODULE  RSCG') 
      IER = 0     
      IF (IPARM(1).LE.0) RETURN       
      N = NN      
      IF (IPARM(11).EQ.0) TIMJ1 = TIMER(DUMMY)  
      IF (LEVEL.GE.3) GO TO 20
      CALL ECHOUT (IPARM,RPARM,6)     
      GO TO 30    
   20 CALL ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,1) 
   30 TEMP = 5.0D2*DRELPR   
      IF (ZETA.GE.TEMP) GO TO 50      
      IF (LEVEL.GE.1) WRITE (NOUT,40) ZETA,DRELPR,TEMP    
   40 FORMAT ('0','*** W A R N I N G ************'/'0',   &
     &   '    IN ITPACK ROUTINE RSCG'/' ','    RPARM(1) =',D10.3,   &
     &   ' (ZETA)'/' ','    A VALUE THIS SMALL MAY HINDER CONVERGENCE '/&
     &   ' ','    SINCE MACHINE PRECISION DRELPR =',D10.3/' ',      &
     &   '    ZETA RESET TO ',D10.3)  
      ZETA = TEMP 
   50 CONTINUE    
      TIME1 = RPARM(9)      
      TIME2 = RPARM(10)     
      DIGIT1 = RPARM(11)    
      DIGIT2 = RPARM(12)    
!       
! ... VERIFY N    
!       
      IF (N.GT.0) GO TO 70  
      IER = 61    
      IF (LEVEL.GE.0) WRITE (NOUT,60) N 
   60 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE RSCG '/' ',      &
     &   '    INVALID MATRIX DIMENSION, N =',I8)
      GO TO 430   
   70 CONTINUE    
!       
! ... REMOVE ROWS AND COLUMNS IF REQUESTED      
!       
      IF (IPARM(10).EQ.0) GO TO 90    
      TOL = RPARM(8)
      CALL IVFILL (N,IWKSP,0) 
      CALL VFILL (N,WKSP,0.0D0)       
      CALL SBELM (N,IA,JA,A,RHS,IWKSP,WKSP,TOL,ISYM,LEVEL,NOUT,IER) 
      IF (IER.EQ.0) GO TO 90
      IF (LEVEL.GE.0) WRITE (NOUT,80) IER,TOL   
   80 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE RSCG '/' ',      &
     &   '    ERROR DETECTED IN SUBROUTINE  SBELM '/' ',  &
     &   '    WHICH REMOVES ROWS AND COLUMNS OF SYSTEM '/' ',       &
     &   '    WHEN DIAGONAL ENTRY TOO LARGE  '/' ','    IER = ',I5,5X,&
     &   ' RPARM(8) = ',D10.3,' (TOL)') 
      GO TO 430   
!       
! ... INITIALIZE WKSP BASE ADDRESSES. 
!       
   90 IB1 = 1     
      IB2 = IB1+N 
      JB3 = IB2+N 
!       
! ... PERMUTE TO  RED-BLACK SYSTEM IF POSSIBLE  
!       
      NB = IPARM(9) 
      IF (NB.GE.0) GO TO 110
      N3 = 3*N    
      CALL IVFILL (N3,IWKSP,0)
      CALL PRBNDX (N,NB,IA,JA,IWKSP,IWKSP(IB2),LEVEL,NOUT,IER)      
      IF (IER.EQ.0) GO TO 110 
      IF (LEVEL.GE.0) WRITE (NOUT,100) IER,NB   
  100 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE RSCG '/' ',      &
     &   '    ERROR DETECTED IN SUBROUTINE  PRBNDX'/' ',  &
     &   '    WHICH COMPUTES THE RED-BLACK INDEXING'/' ','    IER = ',I5&
     &   ,' IPARM(9) = ',I5,' (NB)')  
      GO TO 430   
  110 IF (NB.GE.0.AND.NB.LE.N) GO TO 130
      IER = 64    
      IF (LEVEL.GE.1) WRITE (NOUT,120) IER,NB   
  120 FORMAT (/10X,'ERROR DETECTED IN RED-BLACK SUBSYSTEM INDEX'/10X, &
     &   'IER =',I5,' IPARM(9) =',I5,' (NB)')   
      GO TO 430   
  130 IF (NB.NE.0.AND.NB.NE.N) GO TO 150
      NB = N/2    
      IF (LEVEL.GE.2.AND.IPARM(9).GE.0) WRITE (NOUT,140) IPARM(9),NB
  140 FORMAT (/10X,' IPARM(9) = ',I5,' IMPLIES MATRIX IS DIAGONAL'/10X, &
     &   ' NB RESET TO ',I5)
!       
! ... PERMUTE MATRIX AND RHS
!       
  150 IF (IPARM(9).GE.0) GO TO 190    
      IF (LEVEL.GE.2) WRITE (NOUT,160) NB       
  160 FORMAT (/10X,'ORDER OF BLACK SUBSYSTEM = ',I5,' (NB)')
      CALL PERMAT (N,IA,JA,A,IWKSP,IWKSP(JB3),ISYM,LEVEL,NOUT,IER)  
      IF (IER.EQ.0) GO TO 180 
      IF (LEVEL.GE.0) WRITE (NOUT,170) IER      
  170 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE RSCG '/' ',      &
     &   '    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ',  &
     &   '    WHICH DOES THE RED-BLACK PERMUTATION'/' ','    IER = ',I5)
      GO TO 430   
  180 CALL PERVEC (N,RHS,IWKSP)       
      CALL PERVEC (N,U,IWKSP) 
!       
! ... FINISH WKSP BASE ADDRESSES      
!       
  190 IB3 = IB2+NB
      IB4 = IB3+NB
      IB5 = IB4+NB
      NR = N-NB   
      NRP1 = NR+1 
      IPARM(8) = N+3*NB+2*ITMAX       
      IF (ISYM.NE.0) IPARM(8) = IPARM(8)+2*ITMAX
      IF (NW.GE.IPARM(8)) GO TO 210   
      IER = 62    
      IF (LEVEL.GE.0) WRITE (NOUT,200) NW,IPARM(8)
  200 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE RSCG '/' ',      &
     &   '    NOT ENOUGH WORKSPACE AT ',I10/' ','    SET IPARM(8) =',I10&
     &   ,' (NW)')
      GO TO 430   
!       
! ... SCALE LINEAR SYSTEM, U, AND RHS BY THE SQUARE ROOT OF THE     
! ... DIAGONAL ELEMENTS.    
!       
  210 CONTINUE    
      CALL VFILL (IPARM(8),WKSP,0.0D0)
      CALL SCAL (N,IA,JA,A,RHS,U,WKSP,LEVEL,NOUT,IER)     
      IF (IER.EQ.0) GO TO 230 
      IF (LEVEL.GE.0) WRITE (NOUT,220) IER      
  220 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE RSCG '/' ',      &
     &   '    ERROR DETECTED IN SUBROUTINE  SCAL  '/' ',  &
     &   '    WHICH SCALES THE SYSTEM   '/' ','    IER = ',I5)      
      GO TO 430   
  230 IF (LEVEL.LE.2) GO TO 260       
      WRITE (NOUT,240)      
  240 FORMAT (///1X,'IN THE FOLLOWING, RHO AND GAMMA ARE',&
     &   ' ACCELERATION PARAMETERS')  
      IF (ADAPT) WRITE (NOUT,250)     
  250 FORMAT (1X,'CME IS THE ESTIMATE OF THE LARGEST EIGENVALUE OF',&
     &   ' THE JACOBI MATRIX')
  260 IF (IPARM(11).NE.0) GO TO 270   
      TIMI1 = TIMER(DUMMY)  
!       
! ... INITIALIZE FORWARD PSEUDO-RESIDUAL
!       
  270 CONTINUE    
      IF (N.GT.1) GO TO 280 
      U(1) = RHS(1) 
      GO TO 330   
  280 CALL DCOPY (NR,RHS,1,WKSP(IB1),1) 
      CALL PRSRED (NB,NR,IA,JA,A,U(NRP1),WKSP(IB1))       
      CALL DCOPY (NB,RHS(NRP1),1,WKSP(IB2),1)   
      CALL PRSBLK (NB,NR,IA,JA,A,WKSP(IB1),WKSP(IB2))     
      CALL VEVMW (NB,WKSP(IB2),U(NRP1)) 
!       
! ... ITERATION SEQUENCE    
!       
      ITMAX1 = ITMAX+1      
      DO 300 LOOP = 1,ITMAX1
         IN = LOOP-1
         IF (MOD(IN,2).EQ.1) GO TO 290
!       
! ... CODE FOR THE EVEN ITERATIONS.   
!       
!     U           = U(IN)       WKSP(IB2) = D(IN) 
!     WKSP(IB1)   = U(IN-1)     WKSP(IB3) = D(IN-1)       
!       
         CALL ITRSCG (N,NB,IA,JA,A,U,WKSP(IB1),WKSP(IB2),WKSP(IB3), &
     &      WKSP(IB4),WKSP(IB5))      
!       
         IF (HALT) GO TO 330
         GO TO 300
!       
! ... CODE FOR THE ODD ITERATIONS.    
!       
!     U           = U(IN-1)     WKSP(IB2) = D(IN-1)       
!     WKSP(IB1)   = U(IN)       WKSP(IB3) = D(IN) 
!       
  290    CALL ITRSCG (N,NB,IA,JA,A,WKSP(IB1),U,WKSP(IB3),WKSP(IB2), &
     &      WKSP(IB4),WKSP(IB5))      
!       
         IF (HALT) GO TO 330
  300 CONTINUE    
!       
! ... ITMAX HAS BEEN REACHED
!       
      IF (IPARM(11).NE.0) GO TO 310   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  310 IF (LEVEL.GE.1) WRITE (NOUT,320) ITMAX    
  320 FORMAT ('0','*** W A R N I N G ************'/'0',   &
     &   '    IN ITPACK ROUTINE RSCG'/' ','    FAILURE TO CONVERGE IN', &
     &   I5,' ITERATIONS')  
      IER = 63    
      IF (IPARM(3).EQ.0) RPARM(1) = STPTST      
      GO TO 360   
!       
! ... METHOD HAS CONVERGED  
!       
  330 IF (IPARM(11).NE.0) GO TO 340   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  340 IF (LEVEL.GE.1) WRITE (NOUT,350) IN       
  350 FORMAT (/1X,'RSCG  HAS CONVERGED IN ',I5,' ITERATIONS')       
!       
! ... PUT SOLUTION INTO U IF NOT ALREADY THERE. 
!       
  360 CONTINUE    
      IF (N.EQ.1) GO TO 370 
      IF (MOD(IN,2).EQ.1) CALL DCOPY (N,WKSP(IB1),1,U,1)  
      CALL DCOPY (NR,RHS,1,U,1)       
      CALL PRSRED (NB,NR,IA,JA,A,U(NRP1),U)     
!       
! ... UNSCALE THE MATRIX, SOLUTION, AND RHS VECTORS.      
!       
  370 CALL UNSCAL (N,IA,JA,A,RHS,U,WKSP)
!       
! ... UN-PERMUTE MATRIX,RHS, AND SOLUTION       
!       
      IF (IPARM(9).GE.0) GO TO 400    
      CALL PERMAT (N,IA,JA,A,IWKSP(IB2),IWKSP(JB3),ISYM,LEVEL,NOUT, &
     &   IERPER)  
      IF (IERPER.EQ.0) GO TO 390      
      IF (LEVEL.GE.0) WRITE (NOUT,380) IERPER   
  380 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE RSCG '/' ',      &
     &   '    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ',  &
     &   '    WHICH UNDOES THE RED-BLACK PERMUTATION   '/' ',       &
     &   '    IER = ',I5)   
      IF (IER.EQ.0) IER = IERPER      
      GO TO 430   
  390 CALL PERVEC (N,RHS,IWKSP(IB2))  
      CALL PERVEC (N,U,IWKSP(IB2))    
!       
! ... OPTIONAL ERROR ANALYSIS 
!       
  400 IDGTS = IPARM(12)     
      IF (IDGTS.LT.0) GO TO 410       
      IF (IPARM(2).LE.0) IDGTS = 0    
      CALL PERROR5 (N,IA,JA,A,RHS,U,WKSP,DIGIT1,DIGIT2,IDGTS)
!       
! ... SET RETURN PARAMETERS IN IPARM AND RPARM  
!       
  410 IF (IPARM(11).NE.0) GO TO 420   
      TIMJ2 = TIMER(DUMMY)  
      TIME2 = DBLE(TIMJ2-TIMJ1)       
  420 IPARM(8) = IPARM(8)-2*(ITMAX-IN)
      IF (ISYM.NE.0) IPARM(8) = IPARM(8)-2*(ITMAX-IN)     
      IF (IPARM(3).NE.0) GO TO 430    
      IPARM(1) = IN 
      IPARM(9) = NB 
      RPARM(2) = CME
      RPARM(3) = SME
      RPARM(9) = TIME1      
      RPARM(10) = TIME2     
      RPARM(11) = DIGIT1    
      RPARM(12) = DIGIT2    
!       
  430 CONTINUE    
      IERR = IER  
      IF (LEVEL.GE.3) CALL ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,2)     
!       
      RETURN      
      END 
      SUBROUTINE RSSI (NN,IA,JA,A,RHS,U,IWKSP,NW,WKSP,IPARM,RPARM,IERR) 
!       
!     ITPACK 2C MAIN SUBROUTINE  RSSI  (REDUCED SYSTEM SEMI-ITERATIVE)
!     EACH OF THE MAIN SUBROUTINES:   
!           JCG, JSI, SOR, SSORCG, SSORSI, RSCG, RSSI     
!     CAN BE USED INDEPENDENTLY OF THE OTHERS   
!       
! ... FUNCTION:   
!       
!          THIS SUBROUTINE, RSSI, DRIVES THE  REDUCED SYSTEM SI     
!          ALGORITHM.       
!       
! ... PARAMETER LIST:       
!       
!          N     INPUT INTEGER.  DIMENSION OF THE MATRIX. (= NN)    
!          IA,JA  INPUT INTEGER VECTORS.  THE TWO INTEGER ARRAYS OF 
!                 THE SPARSE MATRIX REPRESENTATION.       
!          A      INPUT D.P. VECTOR.  THE D.P. ARRAY OF THE SPARSE  
!                 MATRIX REPRESENTATION.
!          RHS    INPUT D.P. VECTOR.  CONTAINS THE RIGHT HAND SIDE  
!                 OF THE MATRIX PROBLEM.
!          U      INPUT/OUTPUT D.P. VECTOR.  ON INPUT, U CONTAINS THE 
!                 INITIAL GUESS TO THE SOLUTION. ON OUTPUT, IT CONTAINS 
!                 THE LATEST ESTIMATE TO THE SOLUTION.    
!          IWKSP  INTEGER VECTOR WORKSPACE OF LENGTH 3*N  
!          NW     INPUT INTEGER.  LENGTH OF AVAILABLE WKSP.  ON OUTPUT, 
!                 IPARM(8) IS AMOUNT USED.      
!          WKSP   D.P. VECTOR USED FOR WORKING SPACE.  RSSI 
!                 NEEDS THIS TO BE IN LENGTH AT LEAST  N + NB       
!                 HERE NB IS THE ORDER OF THE BLACK SUBSYSTEM       
!          IPARM  INTEGER VECTOR OF LENGTH 12.  ALLOWS USER TO SPECIFY
!                 SOME INTEGER PARAMETERS WHICH AFFECT THE METHOD.  IF
!          RPARM  D.P. VECTOR OF LENGTH 12. ALLOWS USER TO SPECIFY SOME 
!                 D.P. PARAMETERS WHICH AFFECT THE METHOD.
!          IER     OUTPUT INTEGER.  ERROR FLAG. (= IERR)  
!       
! ... RSSI SUBPROGRAM REFERENCES:     
!       
!          FROM ITPACK    BISRCH, CHEBY, CHGSI, DFAULT, ECHALL,     
!                         ECHOUT, ITERM, TIMER, ITRSSI, IVFILL,     
!                         PARSI, PERMAT, PERROR5, PERVEC, PMULT,     
!                         PRBNDX, PRSBLK, PRSRED, PSTOP, QSORT,     
!                         DAXPY, SBELM, SCAL, DCOPY, DDOT, SUM3,    
!                         TSTCHG, UNSCAL, VEVMW, VFILL, VOUT,       
!                         WEVMW       
!          SYSTEM         DABS, DLOG10, DBLE(AMAX0), DMAX1, DBLE(FLOAT),
!                         DSQRT       
!       
!     VERSION:  ITPACK 2C (MARCH 1982)
!       
!     CODE WRITTEN BY:  DAVID KINCAID, ROGER GRIMES, JOHN RESPESS   
!                       CENTER FOR NUMERICAL ANALYSIS     
!                       UNIVERSITY OF TEXAS     
!                       AUSTIN, TX  78712       
!                       (512) 471-1242
!       
!     FOR ADDITIONAL DETAILS ON THE   
!          (A) SUBROUTINE SEE TOMS ARTICLE 1982 
!          (B) ALGORITHM  SEE CNA REPORT 150    
!       
!     BASED ON THEORY BY:  DAVID YOUNG, DAVID KINCAID, LOU HAGEMAN  
!       
!     REFERENCE THE BOOK:  APPLIED ITERATIVE METHODS      
!                          L. HAGEMAN, D. YOUNG 
!                          ACADEMIC PRESS, 1981 
!       
!     **************************************************  
!     *               IMPORTANT NOTE                   *  
!     *                                                *  
!     *      WHEN INSTALLING ITPACK ROUTINES ON A      *  
!     *  DIFFERENT COMPUTER, RESET SOME OF THE VALUES  *  
!     *  IN  SUBROUTNE DFAULT.   MOST IMPORTANT ARE    *  
!     *                                                *  
!     *   DRELPR      MACHINE RELATIVE PRECISION       *  
!     *   RPARM(1)    STOPPING CRITERION               *  
!     *                                                *  
!     *   ALSO CHANGE SYSTEM-DEPENDENT ROUTINE         *  
!     *   SECOND USED IN TIMER                         *  
!     *                                                *  
!     **************************************************  
!       
!     SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(1),JA(1),IWKSP(1),IPARM(12),NN,NW,IERR   
      DOUBLE PRECISION A(1),RHS(NN),U(NN),WKSP(NW),RPARM(12)
!       
!     SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER IB1,IB2,IDGTS,IER,IERPER,ITMAX1,JB3,LOOP,N,NB,NR,NRP1,N3
      DOUBLE PRECISION DIGIT1,DIGIT2,TEMP,TIME1,TIME2,TOL 
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
! ... VARIABLES IN COMMON BLOCK - ITCOM1
!       
!     IN     - ITERATION NUMBER       
!     IS     - ITERATION NUMBER WHEN PARAMETERS LAST CHANGED
!     ISYM   - SYMMETRIC/NONSYMMETRIC STORAGE FORMAT SWITCH 
!     ITMAX  - MAXIMUM NUMBER OF ITERATIONS ALLOWED       
!     LEVEL  - LEVEL OF OUTPUT CONTROL SWITCH   
!     NOUT   - OUTPUT UNIT NUMBER     
!       
! ... VARIABLES IN COMMON BLOCK - ITCOM2
!       
!     ADAPT  - FULLY ADAPTIVE PROCEDURE SWITCH  
!     BETADT - SWITCH FOR ADAPTIVE DETERMINATION OF BETA  
!     CASEII - ADAPTIVE PROCEDURE CASE SWITCH   
!     HALT   - STOPPING TEST SWITCH   
!     PARTAD - PARTIALLY ADAPTIVE PROCEDURE SWITCH
!       
! ... VARIABLES IN COMMON BLOCK - ITCOM3
!       
!     BDELNM - TWO NORM OF B TIMES DELTA-SUPER-N
!     BETAB  - ESTIMATE FOR THE SPECTRAL RADIUS OF LU MATRIX
!     CME    - ESTIMATE OF LARGEST EIGENVALUE   
!     DELNNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION N      
!     DELSNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION S      
!     FF     - ADAPTIVE PROCEDURE DAMPING FACTOR
!     GAMMA  - ACCELERATION PARAMETER 
!     OMEGA  - OVERRELAXATION PARAMETER FOR SOR AND SSOR  
!     QA     - PSEUDO-RESIDUAL RATIO  
!     QT     - VIRTUAL SPECTRAL RADIUS
!     RHO    - ACCELERATION PARAMETER 
!     RRR    - ADAPTIVE PARAMETER     
!     SIGE   - PARAMETER SIGMA-SUB-E  
!     SME    - ESTIMATE OF SMALLEST EIGENVALUE  
!     SPECR  - SPECTRAL RADIUS ESTIMATE FOR SSOR
!     DRELPR - MACHINE RELATIVE PRECISION       
!     STPTST - STOPPING PARAMETER     
!     UDNM   - TWO NORM OF U
!     ZETA   - STOPPING CRITERION     
!       
! ... INITIALIZE COMMON BLOCKS
!       
      LEVEL = IPARM(2)      
      NOUT = IPARM(4)       
      IF (LEVEL.GE.1) WRITE (NOUT,10) 
   10 FORMAT ('0'///1X,'BEGINNING OF ITPACK SOLUTION MODULE  RSSI') 
      IER = 0     
      IF (IPARM(1).LE.0) RETURN       
      N = NN      
      IF (IPARM(11).EQ.0) TIMJ1 = TIMER(DUMMY)  
      IF (LEVEL.GE.3) GO TO 20
      CALL ECHOUT (IPARM,RPARM,7)     
      GO TO 30    
   20 CALL ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,1) 
   30 TEMP = 5.0D2*DRELPR   
      IF (ZETA.GE.TEMP) GO TO 50      
      IF (LEVEL.GE.1) WRITE (NOUT,40) ZETA,DRELPR,TEMP    
   40 FORMAT ('0','*** W A R N I N G ************'/'0',   &
     &   '    IN ITPACK ROUTINE RSSI'/' ','    RPARM(1) =',D10.3,   &
     &   ' (ZETA)'/' ','    A VALUE THIS SMALL MAY HINDER CONVERGENCE '/&
     &   ' ','    SINCE MACHINE PRECISION DRELPR =',D10.3/' ',      &
     &   '    ZETA RESET TO ',D10.3)  
      ZETA = TEMP 
   50 CONTINUE    
      TIME1 = RPARM(9)      
      TIME2 = RPARM(10)     
      DIGIT1 = RPARM(11)    
      DIGIT2 = RPARM(12)    
!       
! ... VERIFY N    
!       
      IF (N.GT.0) GO TO 70  
      IER = 71    
      IF (LEVEL.GE.0) WRITE (NOUT,60) N 
   60 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE RSSI '/' ',      &
     &   '    INVALID MATRIX DIMENSION, N =',I8)
      GO TO 420   
   70 CONTINUE    
!       
! ... REMOVE ROWS AND COLUMNS IF REQUESTED      
!       
      IF (IPARM(10).EQ.0) GO TO 90    
      TOL = RPARM(8)
      CALL IVFILL (N,IWKSP,0) 
      CALL VFILL (N,WKSP,0.0D0)       
      CALL SBELM (N,IA,JA,A,RHS,IWKSP,WKSP,TOL,ISYM,LEVEL,NOUT,IER) 
      IF (IER.EQ.0) GO TO 90
      IF (LEVEL.GE.0) WRITE (NOUT,80) IER,TOL   
   80 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE RSSI '/' ',      &
     &   '    ERROR DETECTED IN SUBROUTINE  SBELM '/' ',  &
     &   '    WHICH REMOVES ROWS AND COLUMNS OF SYSTEM '/' ',       &
     &   '    WHEN DIAGONAL ENTRY TOO LARGE  '/' ','    IER = ',I5,5X,&
     &   ' RPARM(8) = ',D10.3,' (TOL)') 
!       
! ... INITIALIZE WKSP BASE ADDRESSES. 
!       
   90 IB1 = 1     
      IB2 = IB1+N 
      JB3 = IB2+N 
!       
! ... PERMUTE TO  RED-BLACK SYSTEM IF POSSIBLE  
!       
      NB = IPARM(9) 
      IF (NB.GE.0) GO TO 110
      N3 = 3*N    
      CALL IVFILL (N3,IWKSP,0)
      CALL PRBNDX (N,NB,IA,JA,IWKSP,IWKSP(IB2),LEVEL,NOUT,IER)      
      IF (IER.EQ.0) GO TO 110 
      IF (LEVEL.GE.0) WRITE (NOUT,100) IER,NB   
  100 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE RSSI '/' ',      &
     &   '    ERROR DETECTED IN SUBROUTINE  PRBNDX'/' ',  &
     &   '    WHICH COMPUTES THE RED-BLACK INDEXING'/' ','    IER = ',I5&
     &   ,' IPARM(9) = ',I5,' (NB)')  
      GO TO 420   
  110 IF (NB.GE.0.AND.NB.LE.N) GO TO 130
      IER = 74    
      IF (LEVEL.GE.1) WRITE (NOUT,120) IER,NB   
  120 FORMAT (/10X,'ERROR DETECTED IN RED-BLACK SUBSYSTEM INDEX'/10X, &
     &   'IER =',I5,' IPARM(9) =',I5,' (NB)')   
      GO TO 420   
  130 IF (NB.NE.0.AND.NB.NE.N) GO TO 150
      NB = N/2    
      IF (LEVEL.GE.2.AND.IPARM(9).GE.0) WRITE (NOUT,140) IPARM(9),NB
  140 FORMAT (/10X,' IPARM(9) = ',I5,' IMPLIES MATRIX IS DIAGONAL'/10X, &
     &   ' NB RESET TO ',I5)
!       
! ... PERMUTE MATRIX AND RHS
!       
  150 IF (IPARM(9).GE.0) GO TO 190    
      IF (LEVEL.GE.2) WRITE (NOUT,160) NB       
  160 FORMAT (/10X,'ORDER OF BLACK SUBSYSTEM = ',I5,' (NB)')
      CALL PERMAT (N,IA,JA,A,IWKSP,IWKSP(JB3),ISYM,LEVEL,NOUT,IER)  
      IF (IER.EQ.0) GO TO 180 
      IF (LEVEL.GE.0) WRITE (NOUT,170) IER      
  170 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE RSSI '/' ',      &
     &   '    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ',  &
     &   '    WHICH DOES THE RED-BLACK PERMUTATION'/' ','    IER = ',I5)
      GO TO 420   
  180 CALL PERVEC (N,RHS,IWKSP)       
      CALL PERVEC (N,U,IWKSP) 
!       
! ... INITIALIZE WKSP BASE ADDRESSES  
!       
  190 NR = N-NB   
!       
      NRP1 = NR+1 
      IPARM(8) = N+NB       
      IF (NW.GE.IPARM(8)) GO TO 210   
      IER = 72    
      IF (LEVEL.GE.0) WRITE (NOUT,200) NW,IPARM(8)
  200 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE RSSI '/' ',      &
     &   '    NOT ENOUGH WORKSPACE AT ',I10/' ','    SET IPARM(8) =',I10&
     &   ,' (NW)')
      GO TO 420   
!       
! ... SCALE LINEAR SYSTEM, U, AND RHS BY THE SQUARE ROOT OF THE     
! ... DIAGONAL ELEMENTS.    
!       
  210 CONTINUE    
      CALL VFILL (IPARM(8),WKSP,0.0D0)
      CALL SCAL (N,IA,JA,A,RHS,U,WKSP,LEVEL,NOUT,IER)     
      IF (IER.EQ.0) GO TO 230 
      IF (LEVEL.GE.0) WRITE (NOUT,220) IER      
  220 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE RSSI '/' ',      &
     &   '    ERROR DETECTED IN SUBROUTINE  SCAL  '/' ',  &
     &   '    WHICH SCALES THE SYSTEM   '/' ','    IER = ',I5)      
      GO TO 420   
  230 IF (LEVEL.LE.2) GO TO 250       
      WRITE (NOUT,240)      
  240 FORMAT (///1X,'IN THE FOLLOWING, RHO AND GAMMA ARE',&
     &   ' ACCELERATION PARAMETERS')  
  250 IF (IPARM(11).NE.0) GO TO 260   
      TIMI1 = TIMER(DUMMY)  
!       
! ... ITERATION SEQUENCE    
!       
  260 IF (N.GT.1) GO TO 270 
      U(1) = RHS(1) 
      GO TO 320   
  270 ITMAX1 = ITMAX+1      
      DO 290 LOOP = 1,ITMAX1
         IN = LOOP-1
         IF (MOD(IN,2).EQ.1) GO TO 280
!       
! ... CODE FOR THE EVEN ITERATIONS.   
!       
!     U           = U(IN)   
!     WKSP(IB1)   = U(IN-1) 
!       
         CALL ITRSSI (N,NB,IA,JA,A,RHS,U,WKSP(IB1),WKSP(IB2))       
!       
         IF (HALT) GO TO 320
         GO TO 290
!       
! ... CODE FOR THE ODD ITERATIONS.    
!       
!     U           = U(IN-1) 
!     WKSP(IB1)   = U(IN)   
!       
  280    CALL ITRSSI (N,NB,IA,JA,A,RHS,WKSP(IB1),U,WKSP(IB2))       
!       
         IF (HALT) GO TO 320
  290 CONTINUE    
!       
! ... ITMAX HAS BEEN REACHED
!       
      IF (IPARM(11).NE.0) GO TO 300   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  300 IF (LEVEL.GE.1) WRITE (NOUT,310) ITMAX    
  310 FORMAT ('0','*** W A R N I N G ************'/'0',   &
     &   '    IN ITPACK ROUTINE RSSI'/' ','    FAILURE TO CONVERGE IN', &
     &   I5,' ITERATIONS')  
      IER = 73    
      IF (IPARM(3).EQ.0) RPARM(1) = STPTST      
      GO TO 350   
!       
! ... METHOD HAS CONVERGED  
!       
  320 IF (IPARM(11).NE.0) GO TO 330   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  330 IF (LEVEL.GE.1) WRITE (NOUT,340) IN       
  340 FORMAT (/1X,'RSSI  HAS CONVERGED IN ',I5,' ITERATIONS')       
!       
! ... PUT SOLUTION INTO U IF NOT ALREADY THERE. 
!       
  350 CONTINUE    
      IF (N.EQ.1) GO TO 360 
      IF (MOD(IN,2).EQ.1) CALL DCOPY (N,WKSP(IB1),1,U,1)  
      CALL DCOPY (NR,RHS,1,U,1)       
      CALL PRSRED (NB,NR,IA,JA,A,U(NRP1),U)     
!       
! ... UNSCALE THE MATRIX, SOLUTION, AND RHS VECTORS.      
!       
  360 CALL UNSCAL (N,IA,JA,A,RHS,U,WKSP)
!       
! ... UN-PERMUTE MATRIX,RHS, AND SOLUTION       
!       
      IF (IPARM(9).GE.0) GO TO 390    
      CALL PERMAT (N,IA,JA,A,IWKSP(IB2),IWKSP(JB3),ISYM,LEVEL,NOUT, &
     &   IERPER)  
      IF (IERPER.EQ.0) GO TO 380      
      IF (LEVEL.GE.0) WRITE (NOUT,370) IERPER   
  370 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    CALLED FROM ITPACK ROUTINE RSSI '/' ',      &
     &   '    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ',  &
     &   '    WHICH UNDOES THE RED-BLACK PERMUTATION   '/' ',       &
     &   '    IER = ',I5)   
      IF (IER.EQ.0) IER = IERPER      
      GO TO 420   
  380 CALL PERVEC (N,RHS,IWKSP(IB2))  
      CALL PERVEC (N,U,IWKSP(IB2))    
!       
! ... OPTIONAL ERROR ANALYSIS 
!       
  390 IDGTS = IPARM(12)     
      IF (IDGTS.LT.0) GO TO 400       
      IF (IPARM(2).LE.0) IDGTS = 0    
      CALL PERROR5 (N,IA,JA,A,RHS,U,WKSP,DIGIT1,DIGIT2,IDGTS)
!       
! ... SET RETURN PARAMETERS IN IPARM AND RPARM  
!       
  400 IF (IPARM(11).NE.0) GO TO 410   
      TIMJ2 = TIMER(DUMMY)  
      TIME2 = DBLE(TIMJ2-TIMJ1)       
  410 IF (IPARM(3).NE.0) GO TO 420    
      IPARM(1) = IN 
      IPARM(9) = NB 
      RPARM(2) = CME
      RPARM(3) = SME
      RPARM(9) = TIME1      
      RPARM(10) = TIME2     
      RPARM(11) = DIGIT1    
      RPARM(12) = DIGIT2    
!       
  420 CONTINUE    
      IERR = IER  
      IF (LEVEL.GE.3) CALL ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,2)     
!       
      RETURN      
      END













      SUBROUTINE ITJCG (NN,IA,JA,A,U,U1,D,D1,DTWD,TRI)    
!       
! ... FUNCTION:   
!       
!          THIS SUBROUTINE, ITJCG, PERFORMS ONE ITERATION OF THE    
!          JACOBI CONJUGATE GRADIENT ALGORITHM.  IT IS CALLED BY JCG. 
!       
! ... PARAMETER LIST:       
!       
!          N      INPUT INTEGER.  DIMENSION OF THE MATRIX. (= NN)   
!          IA,JA  INPUT INTEGER VECTORS.  CONTAINS INFORMATION DEFINING 
!                 THE SPARSE MATRIX REPRESENTATION.       
!          A      INPUT D.P. VECTOR. CONTAINS THE NONZERO VALUES OF THE 
!                 LINEAR SYSTEM.      
!          U      INPUT D.P. VECTOR.  CONTAINS THE VALUE OF THE     
!                 SOLUTION VECTOR AT THE END OF IN ITERATIONS.      
!          U1     INPUT/OUTPUT D.P. VECTOR.  ON INPUT, IT CONTAINS  
!                 THE VALUE OF THE SOLUTION AT THE END OF THE IN-1  
!                 ITERATION.  ON OUTPUT, IT WILL CONTAIN THE NEWEST 
!                 ESTIMATE FOR THE SOLUTION VECTOR.       
!          D      INPUT D.P. VECTOR.  CONTAINS THE PSEUDO-RESIDUAL  
!                 VECTOR AFTER IN ITERATIONS.   
!          D1     INPUT/OUTPUT D.P. VECTOR.  ON INPUT, D1 CONTAINS  
!                 THE PSEUDO-RESIDUAL VECTOR AFTER IN-1 ITERATIONS.  ON 
!                 OUTPUT, IT WILL CONTAIN THE NEWEST PSEUDO-RESIDUAL
!                 VECTOR.   
!          DTWD   D.P. ARRAY.  USED IN THE COMPUTATIONS OF THE      
!                 ACCELERATION PARAMETER GAMMA AND THE NEW PSEUDO-  
!                 RESIDUAL. 
!          TRI    D.P. ARRAY.  STORES THE TRIDIAGONAL MATRIX ASSOCIATED 
!                 WITH THE EIGENVALUES OF THE CONJUGATE GRADIENT    
!                 POLYNOMIAL. 
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(1),JA(1),NN
      DOUBLE PRECISION A(1),U(NN),U1(NN),D(NN),D1(NN),DTWD(NN),TRI(2,1) 
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER N   
      DOUBLE PRECISION CON,C1,C2,C3,C4,DNRM,DTNRM,GAMOLD,RHOOLD,RHOTMP
      LOGICAL Q1  
!       
! ... SPECIFICATIONS FOR FUNCTION SUBPROGRAMS   
!       
      DOUBLE PRECISION DDOT 
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN SUBROUTINE JCG   
!       
! ... COMPUTE NEW ESTIMATE FOR CME IF ADAPT = .TRUE.
!
      IF (ADAPT) CALL CHGCON (TRI,GAMOLD,RHOOLD,1)
!       
! ... TEST FOR STOPPING     
!       
      N = NN      
      DELNNM = DDOT(N,D,1,D,1)
      DNRM = DELNNM 
      CON = CME   
      CALL PSTOP (N,U,DNRM,CON,1,Q1)  
      IF (HALT) GO TO 30    
!       
! ... COMPUTE RHO AND GAMMA - ACCELERATION PARAMETERS     
!       
      CALL VFILL (N,DTWD,0.D0)
      CALL PJAC (N,IA,JA,A,D,DTWD)    
      DTNRM = DDOT(N,D,1,DTWD,1)      
      IF (ISYM.EQ.0) GO TO 10 
      RHOTMP = DDOT(N,DTWD,1,D1,1)    
      CALL PARCON (DTNRM,C1,C2,C3,C4,GAMOLD,RHOTMP,1)     
      RHOOLD = RHOTMP       
      GO TO 20    
   10 CALL PARCON (DTNRM,C1,C2,C3,C4,GAMOLD,RHOOLD,1)     
!       
! ... COMPUTE U(IN+1) AND D(IN+1)     
!       
   20 CALL SUM3 (N,C1,D,C2,U,C3,U1)   
      CALL SUM3 (N,C1,DTWD,C4,D,C3,D1)
!       
! ... OUTPUT INTERMEDIATE INFORMATION 
!       
   30 CALL ITERM (N,A,U,DTWD,1)       
!       
      RETURN      
      END















      SUBROUTINE ITJSI (NN,IA,JA,A,RHS,U,U1,D,ICNT)       
!       
! ... FUNCTION:   
!       
!          THIS SUBROUTINE, ITJSI, PERFORMS ONE ITERATION OF THE    
!          JACOBI SEMI-ITERATIVE ALGORITHM.  IT IS CALLED BY JSI.   
!       
! ... PARAMETER LIST:       
!       
!          N      INPUT INTEGER.  DIMENSION OF THE MATRIX. (= NN)   
!          IA,JA  INPUT INTEGER VECTORS.  THE TWO INTEGER ARRAYS OF 
!                 THE SPARSE MATRIX REPRESENTATION.       
!          A      INPUT D.P. VECTOR.  THE D.P. ARRAY OF THE SPARSE  
!                 MATRIX REPRESENTATION.
!          RHS    INPUT D.P. VECTOR.  CONTAINS THE RIGHT HAND SIDE  
!                 OF THE MATRIX PROBLEM.
!          U      INPUT D.P. VECTOR.  CONTAINS THE ESTIMATE FOR THE 
!                 SOLUTION VECTOR AFTER IN ITERATIONS.    
!          U1     INPUT/OUTPUT D.P. VECTOR.  ON INPUT, U1 CONTAINS THE
!                 SOLUTION VECTOR AFTER IN-1 ITERATIONS.  ON OUTPUT,
!                 IT WILL CONTAIN THE NEWEST ESTIMATE FOR THE SOLUTION
!                 VECTOR.   
!          D      D.P. ARRAY.  D IS USED FOR THE COMPUTATION OF THE 
!                 PSEUDO-RESIDUAL ARRAY FOR THE CURRENT ITERATION.  
!          ICNT   NUMBER OF ITERATIONS SINCE LAST CHANGE OF SME     
!       
! ... SPECIFICATIONS OF ARGUMENTS     
!       
      INTEGER IA(1),JA(1),NN,ICNT     
      DOUBLE PRECISION A(1),RHS(NN),U(NN),U1(NN),D(NN)    
!       
! ... SPECIFICATIONS OF LOCAL VARIABLES 
!       
      INTEGER N   
      DOUBLE PRECISION CON,C1,C2,C3,DNRM,DTNRM,OLDNRM     
      LOGICAL Q1  
!       
! ... SPECIFICATIONS OF FUNCTION SUBPROGRAMS    
!       
      DOUBLE PRECISION DDOT,PVTBV     
      LOGICAL TSTCHG,CHGSME 
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN SUBROUTINE JSI   
!       
      N = NN      
      IF (IN.EQ.0) ICNT = 0 
!       
! ... COMPUTE PSEUDO-RESIDUALS
!       
      CALL DCOPY (N,RHS,1,D,1)
      CALL PJAC (N,IA,JA,A,U,D)       
      CALL VEVMW (N,D,U)    
!       
! ... STOPPING AND ADAPTIVE CHANGE TESTS
!       
      OLDNRM = DELNNM       
      DELNNM = DDOT(N,D,1,D,1)
      DNRM = DELNNM 
      CON = CME   
      CALL PSTOP (N,U,DNRM,CON,1,Q1)  
      IF (HALT) GO TO 40    
      IF (.NOT.ADAPT) GO TO 30
      IF (.NOT.TSTCHG(1)) GO TO 10    
!       
! ... CHANGE ITERATIVE PARAMETERS (CME) 
!       
      DTNRM = PVTBV(N,IA,JA,A,D)      
      CALL CHGSI (DTNRM,1)  
      IF (.NOT.ADAPT) GO TO 30
      GO TO 20    
!       
! ... TEST IF SME NEEDS TO BE CHANGED AND CHANGE IF NECESSARY.      
!       
   10 CONTINUE    
      IF (CASEII) GO TO 30  
      IF (.NOT.CHGSME(OLDNRM,ICNT)) GO TO 30    
      ICNT = 0    
!       
! ... COMPUTE U(IN+1) AFTER CHANGE OF PARAMETERS
!       
   20 CALL DCOPY (N,U,1,U1,1) 
      CALL DAXPY (N,GAMMA,D,1,U1,1)   
      GO TO 40    
!       
! ... COMPUTE U(IN+1) WITHOUT CHANGE OF PARAMETERS
!       
   30 CALL PARSI (C1,C2,C3,1) 
      CALL SUM3 (N,C1,D,C2,U,C3,U1)   
!       
! ... OUTPUT INTERMEDIATE INFORMATION 
!       
   40 CALL ITERM (N,A,U,D,2)
!       
      RETURN      
      END 
      SUBROUTINE ITSOR (NN,IA,JA,A,RHS,U,WK)    
!       
! ... FUNCTION:   
!       
!          THIS SUBROUTINE, ITSOR, PERFORMS ONE ITERATION OF THE    
!          SUCCESSIVE OVERRELAXATION ALGORITHM.  IT IS CALLED BY SOR. 
!       
! ... PARAMETER LIST:       
!       
!          N      INPUT INTEGER.  DIMENSION OF THE MATRIX. (= NN)   
!          IA,JA  INPUT INTEGER VECTORS.  THE TWO INTEGER ARRAYS OF 
!                 THE SPARSE MATRIX REPRESENTATION.       
!          A      INPUT D.P. VECTOR.  THE D.P. ARRAY OF THE SPARSE  
!                 MATRIX REPRESENTATION.
!          RHS    INPUT D.P. VECTOR.  CONTAINS THE RIGHT HAND SIDE  
!                 OF THE MATRIX PROBLEM.
!          U      INPUT/OUTPUT D.P. VECTOR.  ON INPUT, U CONTAINS THE 
!                 SOLUTION VECTOR AFTER IN ITERATIONS.  ON OUTPUT,  
!                 IT WILL CONTAIN THE NEWEST ESTIMATE FOR THE SOLUTION
!                 VECTOR.   
!          WK     D.P. ARRAY.  WORK VECTOR OF LENGTH N.   
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(1),JA(1),NN
      DOUBLE PRECISION A(1),RHS(NN),U(NN),WK(NN)
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER IP,IPHAT,IPSTAR,ISS,N   
      DOUBLE PRECISION DNRM,H,OMEGAP,SPCRM1     
      LOGICAL CHANGE,Q1     
!       
      DOUBLE PRECISION TAU  
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN SUBROUTINE SOR   
!       
! ... SET INITIAL PARAMETERS NOT ALREADY SET    
!       
      N = NN      
      IF (IN.NE.0) GO TO 20 
      CALL PSTOP (N,U,0.D0,0.D0,0,Q1) 
      IF (ADAPT) GO TO 10   
      CHANGE = .FALSE.      
      IP = 0      
      IPHAT = 2   
      ISS = 0     
      GO TO 30    
!       
   10 CHANGE = .TRUE.       
      IP = 0      
      OMEGAP = OMEGA
      OMEGA = 1.D0
      ISS = 0     
      IPHAT = 2   
      IPSTAR = 4  
      IF (OMEGAP.LE.1.D0) CHANGE = .FALSE.      
!       
! ... RESET OMEGA, IPHAT, AND IPSTAR (CIRCLE A IN FLOWCHART)
!       
   20 IF (.NOT.CHANGE) GO TO 30       
      CHANGE = .FALSE.      
      IS = IS+1   
      IP = 0      
      ISS = 0     
      OMEGA = DMIN1(OMEGAP,TAU(IS))   
      IPHAT = MAX0(3,IFIX(SNGL((OMEGA-1.D0)/(2.D0-OMEGA)))) 
      IPSTAR = IPSTR(OMEGA) 
!       
! ... COMPUTE U (IN + 1) AND NORM OF DEL(S,P) - CIRCLE B IN FLOW CHART
!       
   30 CONTINUE    
      DELSNM = DELNNM       
      SPCRM1 = SPECR
      CALL DCOPY (N,RHS,1,WK,1)       
      CALL PFSOR1 (N,IA,JA,A,U,WK)    
      IF (DELNNM.EQ.0.D0) GO TO 40    
      IF (IN.NE.0) SPECR = DELNNM/DELSNM
      IF (IP.LT.IPHAT) GO TO 70       
!       
! ... STOPPING TEST, SET H  
!       
      IF (SPECR.GE.1.D0) GO TO 70     
      IF (.NOT.(SPECR.GT.(OMEGA-1.D0))) GO TO 40
      H = SPECR   
      GO TO 50    
   40 ISS = ISS+1 
      H = OMEGA-1.D0
!       
! ... PERFORM STOPPING TEST.
!       
   50 CONTINUE    
      DNRM = DELNNM**2      
      CALL PSTOP (N,U,DNRM,H,1,Q1)    
      IF (HALT) GO TO 70    
!       
! ... METHOD HAS NOT CONVERGED YET, TEST FOR CHANGING OMEGA 
!       
      IF (.NOT.ADAPT) GO TO 70
      IF (IP.LT.IPSTAR) GO TO 70      
      IF (OMEGA.GT.1.D0) GO TO 60     
      CME = DSQRT(DABS(SPECR))
      OMEGAP = 2.D0/(1.D0+DSQRT(DABS(1.D0-SPECR)))
      CHANGE = .TRUE.       
      GO TO 70    
   60 IF (ISS.NE.0) GO TO 70
      IF (SPECR.LE.(OMEGA-1.D0)**FF) GO TO 70   
      IF ((SPECR+5.D-5).LE.SPCRM1) GO TO 70     
!       
! ... CHANGE PARAMETERS     
!       
      CME = (SPECR+OMEGA-1.D0)/(DSQRT(DABS(SPECR))*OMEGA) 
      OMEGAP = 2.D0/(1.D0+DSQRT(DABS(1.D0-CME*CME)))      
      CHANGE = .TRUE.       
!       
! ... OUTPUT INTERMEDIATE INFORMATION 
!       
   70 CALL ITERM (N,A,U,WK,3) 
      IP = IP+1   
!       
      RETURN      
      END 
      SUBROUTINE ITSRCG (NN,IA,JA,A,RHS,U,U1,C,C1,D,DL,WK,TRI)      
!       
! ... FUNCTION:   
!       
!          THIS SUBROUTINE, ITSRCG, PERFORMS ONE ITERATION OF THE   
!          SYMMETRIC SOR CONJUGATE GRADIENT ALGORITHM.  IT IS CALLED BY 
!          SSORCG.
!       
! ... PARAMETER LIST:       
!       
!          N      INPUT INTEGER.  DIMENSION OF THE MATRIX. (= NN)   
!          IA,JA  INPUT INTEGER VECTORS.  THE TWO INTEGER ARRAYS OF 
!                 THE SPARSE MATRIX REPRESENTATION.       
!          A      INPUT D.P. VECTOR.  THE D.P. ARRAY OF THE SPARSE  
!                 MATRIX REPRESENTATION.
!          RHS    INPUT D.P. VECTOR.  CONTAINS THE RIGHT HAND SIDE  
!                 OF THE MATRIX PROBLEM.
!          U      INPUT D.P. VECTOR.  CONTAINS THE ESTIMATE OF THE  
!                 SOLUTION VECTOR AFTER IN ITERATIONS.    
!          U1     INPUT/OUTPUT D.P. VECTOR.  ON INPUT, U1 CONTAINS THE
!                 THE ESTIMATE FOR THE SOLUTION AFTER IN-1 ITERATIONS.
!                 ON OUTPUT, U1 CONTAINS THE UPDATED ESTIMATE.      
!          C      INPUT D.P. VECTOR.  CONTAINS THE FORWARD RESIDUAL 
!                 AFTER IN ITERATIONS.
!          C1     INPUT/OUTPUT D.P. VECTOR.  ON INPUT, C1 CONTAINS  
!                 THE FORWARD RESIDUAL AFTER IN-1 ITERATIONS.  ON   
!                 OUTPUT, C1 CONTAINS THE UPDATED FORWARD RESIDUAL. 
!          D      D.P. VECTOR.  IS USED TO COMPUTE THE BACKWARD PSEUDO- 
!                 RESIDUAL VECTOR FOR THE CURRENT ITERATION.
!          DL     D.P. VECTOR.  IS USED IN THE COMPUTATIONS OF THE  
!                 ACCELERATION PARAMETERS.      
!          WK     D.P. VECTOR.  WORKING SPACE OF LENGTH N.
!          TRI    D.P. VECTOR. STORES THE TRIDIAGONAL MATRIX ASSOCIATED 
!                 WITH THE CONJUGATE GRADIENT ACCELERATION. 
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(1),JA(1),NN
      DOUBLE PRECISION A(1),RHS(NN),U(NN),U1(NN),C(NN),C1(NN),D(NN),&
     &   DL(NN),WK(NN),TRI(2,1)       
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER N   
      DOUBLE PRECISION BETNEW,CON,DNRM,GAMOLD,RHOOLD,RHOTMP,T1,T2,T3,T4 
      LOGICAL Q1  
!       
! ... SPECIFICATIONS FOR FUNCTION SUBPROGRAMS   
!       
      DOUBLE PRECISION DDOT,PBETA,PVTBV 
      LOGICAL OMGCHG,OMGSTR 
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN SUBROUTINE SSORCG
!       
! ... CALCULATE S-PRIME FOR ADAPTIVE PROCEDURE. 
!       
      N = NN      
      IF (ADAPT.OR.PARTAD) CALL CHGCON (TRI,GAMOLD,RHOOLD,3)
!       
! ... COMPUTE BACKWARD RESIDUAL       
!       
      CALL DCOPY (N,RHS,1,WK,1)       
      CALL DCOPY (N,C,1,D,1)
      CALL VEVPW (N,D,U)    
      CALL PBSOR (N,IA,JA,A,D,WK)     
      CALL VEVMW (N,D,U)    
!       
! ... COMPUTE ACCELERATION PARAMETERS AND THEN U(IN+1) (IN U1)      
!       
      CALL DCOPY (N,D,1,DL,1) 
      CALL VFILL (N,WK,0.D0)
      CALL PFSOR (N,IA,JA,A,DL,WK)    
      CALL WEVMW (N,D,DL)   
      DELNNM = DDOT(N,C,1,C,1)
      IF (DELNNM.EQ.0.D0) GO TO 30    
      DNRM = DDOT(N,C,1,DL,1) 
      IF (DNRM.EQ.0.D0) GO TO 30      
      IF (ISYM.EQ.0) GO TO 10 
      RHOTMP = DDOT(N,C,1,C1,1)-DDOT(N,DL,1,C1,1) 
      CALL PARCON (DNRM,T1,T2,T3,T4,GAMOLD,RHOTMP,3)      
      RHOOLD = RHOTMP       
      GO TO 20    
   10 CALL PARCON (DNRM,T1,T2,T3,T4,GAMOLD,RHOOLD,3)      
   20 CALL SUM3 (N,T1,D,T2,U,T3,U1)   
!       
! ... TEST FOR STOPPING     
!       
   30 BDELNM = DDOT(N,D,1,D,1)
      DNRM = BDELNM 
      CON = SPECR 
      CALL PSTOP (N,U,DNRM,CON,1,Q1)  
      IF (HALT) GO TO 100   
!       
! ... IF NON- OR PARTIALLY-ADAPTIVE, COMPUTE C(IN+1) AND EXIT.      
!       
      IF (ADAPT) GO TO 40   
      CALL SUM3 (N,-T1,DL,T2,C,T3,C1) 
      GO TO 100   
!       
! ... FULLY ADAPTIVE PROCEDURE
!       
   40 CONTINUE    
      IF (OMGSTR(1)) GO TO 90 
      IF (OMGCHG(1)) GO TO 50 
!       
! ... PARAMETERS HAVE BEEN UNCHANGED.  COMPUTE C(IN+1) AND EXIT.    
!       
      CALL SUM3 (N,-T1,DL,T2,C,T3,C1) 
      GO TO 100   
!       
! ... IT HAS BEEN DECIDED TO CHANGE PARAMETERS  
!        (1) COMPUTE NEW BETAB IF BETADT = .TRUE. 
!       
   50 CONTINUE    
      IF (.NOT.BETADT) GO TO 60       
      BETNEW = PBETA(N,IA,JA,A,D,WK,C1)/BDELNM  
      BETAB = DMAX1(BETAB,.25D0,BETNEW) 
!       
! ...    (2) COMPUTE NEW CME, OMEGA, AND SPECR  
!       
   60 CONTINUE    
      IF (CASEII) GO TO 70  
      DNRM = PVTBV(N,IA,JA,A,D)       
      GO TO 80    
   70 CALL VFILL (N,WK,0.D0)
      CALL PJAC (N,IA,JA,A,D,WK)      
      DNRM = DDOT(N,WK,1,WK,1)
   80 CALL OMEG (DNRM,3)    
!       
! ...    (3) COMPUTE NEW FORWARD RESIDUAL SINCE OMEGA HAS BEEN CHANGED. 
!       
   90 CONTINUE    
      CALL DCOPY (N,RHS,1,WK,1)       
      CALL DCOPY (N,U1,1,C1,1)
      CALL PFSOR (N,IA,JA,A,C1,WK)    
      CALL VEVMW (N,C1,U1)  
!       
! ... OUTPUT INTERMEDIATE RESULTS.    
!       
  100 CALL ITERM (N,A,U,WK,4) 
!       
      RETURN      
      END 
      SUBROUTINE ITSRSI (NN,IA,JA,A,RHS,U,U1,C,D,CTWD,WK) 
!       
! ... FUNCTION:   
!       
!          THIS SUBROUTINE, ITSRSI, PERFORMS ONE ITERATION OF THE   
!          SYMMETRIC SOR SEMI-ITERATION ALGORITHM.  IT IS CALLED BY 
!          SSORSI.
!       
! ... PARAMETER LIST:       
!       
!          N      INPUT INTEGER.  DIMENSION OF THE MATRIX. (= NN)   
!          IA,JA  INPUT INTEGER VECTORS.  THE TWO INTEGER ARRAYS OF 
!                 THE SPARSE MATRIX REPRESENTATION.       
!          A      INPUT D.P. VECTOR.  THE D.P. ARRAY OF THE SPARSE  
!                 MATRIX REPRESENTATION.
!          RHS    INPUT D.P. VECTOR.  CONTAINS THE RIGHT HAND SIDE  
!                 OF THE MATRIX PROBLEM.
!          U      INPUT D.P. VECTOR.  CONTAINS THE ESTIMATE OF THE  
!                 SOLUTION VECTOR AFTER IN ITERATIONS.    
!          U1     INPUT/OUTPUT D.P. VECTOR.  ON INPUT, U1 CONTAINS THE
!                 THE ESTIMATE FOR THE SOLUTION AFTER IN-1 ITERATIONS.
!                 ON OUTPUT, U1 CONTAINS THE UPDATED ESTIMATE.      
!          C      D.P. VECTOR.  IS USED TO COMPUTE THE FORWARD PSEUDO-
!                 RESIDUAL VECTOR FOR THE CURRENT ITERATION.
!          D      D.P. VECTOR.  IS USED TO COMPUTE THE BACKWARD PSEUDO- 
!                 RESIDUAL VECTOR FOR THE CURRENT ITERATION.
!          CTWD   D.P. VECTOR.  IS USED IN THE COMPUTATIONS OF THE  
!                 ACCELERATION PARAMETERS.      
!          WK     D.P. VECTOR.  WORKING SPACE OF LENGTH N.
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(1),JA(1),NN
      DOUBLE PRECISION A(1),RHS(NN),U(NN),U1(NN),C(NN),D(NN),WK(NN),&
     &   CTWD(NN) 
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER N   
      DOUBLE PRECISION BETNEW,CON,C1,C2,C3,DNRM 
      LOGICAL Q1  
!       
! ... SPECIFICATIONS FOR FUNCTION SUBPROGRAMS   
!       
      DOUBLE PRECISION DDOT,PBETA,PVTBV 
      LOGICAL OMGSTR,TSTCHG 
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN SUBROUTINE SSORSI
!       
! ... COMPUTE PSEUDO-RESIDUALS (FORWARD AND BACKWARD)     
!       
      N = NN      
      CALL DCOPY (N,RHS,1,WK,1)       
      CALL DCOPY (N,U,1,CTWD,1)       
      CALL PSSOR1 (N,IA,JA,A,CTWD,WK,C,D)       
!       
! ... COMPUTE U(IN+1) -- CONTAINED IN THE VECTOR U1.      
!       
      CALL PARSI (C1,C2,C3,3) 
      CALL SUM3 (N,C1,D,C2,U,C3,U1)   
!       
! ... TEST FOR STOPPING     
!       
      BDELNM = DDOT(N,D,1,D,1)
      DNRM = BDELNM 
      CON = SPECR 
      CALL PSTOP (N,U,DNRM,CON,1,Q1)  
      IF (HALT.OR..NOT.(ADAPT.OR.PARTAD)) GO TO 40
!       
! ... ADAPTIVE PROCEDURE    
!       
      IF (OMGSTR(1)) GO TO 40 
      DELNNM = DDOT(N,C,1,C,1)
      IF (IN.EQ.IS) DELSNM = DELNNM   
      IF (IN.EQ.0.OR..NOT.TSTCHG(1)) GO TO 40   
!       
! ... IT HAS BEEN DECIDED TO CHANGE PARAMETERS. 
! ...    (1) COMPUTE CTWD   
!       
      CALL DCOPY (N,D,1,CTWD,1)       
      CALL VFILL (N,WK,0.D0)
      CALL PFSOR (N,IA,JA,A,CTWD,WK)  
      CALL VEVPW (N,CTWD,C) 
      CALL VEVMW (N,CTWD,D) 
!       
! ...    (2) COMPUTE NEW SPECTRAL RADIUS FOR CURRENT OMEGA. 
!       
      DNRM = DDOT(N,C,1,CTWD,1)       
      CALL CHGSI (DNRM,3)   
      IF (.NOT.ADAPT) GO TO 40
!       
! ...    (3) COMPUTE NEW BETAB IF BETADT = .TRUE. 
!       
      IF (.NOT.BETADT) GO TO 10       
      BETNEW = PBETA(N,IA,JA,A,D,WK,CTWD)/BDELNM
      BETAB = DMAX1(BETAB,.25D0,BETNEW) 
!       
! ...    (4) COMPUTE NEW CME, OMEGA, AND SPECR. 
!       
   10 CONTINUE    
      IF (CASEII) GO TO 20  
      DNRM = PVTBV(N,IA,JA,A,D)       
      GO TO 30    
   20 CALL VFILL (N,WK,0.D0)
      CALL PJAC (N,IA,JA,A,D,WK)      
      DNRM = DDOT(N,WK,1,WK,1)
   30 CALL OMEG (DNRM,3)    
!       
! ... OUTPUT INTERMEDIATE INFORMATION 
!       
   40 CALL ITERM (N,A,U,WK,5) 
!       
      RETURN      
      END 
      SUBROUTINE ITRSCG (N,NNB,IA,JA,A,UB,UB1,DB,DB1,WB,TRI)
!       
! ... FUNCTION:   
!       
!          THIS SUBROUTINE, ITRSCG, PERFORMS ONE ITERATION OF THE   
!          REDUCED SYSTEM CONJUGATE GRADIENT ALGORITHM.  IT IS      
!          CALLED BY RSCG.  
!       
! ... PARAMETER LIST:       
!       
!          N      INPUT INTEGER.  DIMENSION OF THE MATRIX.
!          NB     INPUT INTEGER.  CONTAINS THE NUMBER OF BLACK POINTS 
!                 IN THE RED-BLACK MATRIX. (= NNB)
!          IA,JA  INPUT INTEGER VECTORS.  THE TWO INTEGER ARRAYS OF 
!                 THE SPARSE MATRIX REPRESENTATION.       
!          A      INPUT D.P. VECTOR.  THE D.P. ARRAY OF THE SPARSE  
!                 MATRIX REPRESENTATION.
!          UB     INPUT D.P. VECTOR.  CONTAINS THE ESTIMATE FOR THE 
!                 SOLUTION ON THE BLACK POINTS AFTER IN ITERATIONS. 
!          UB1    INPUT/OUTPUT D.P. VECTOR.  ON INPUT, UB1 CONTAINS THE 
!                 SOLUTION VECTOR AFTER IN-1 ITERATIONS.  ON OUTPUT,
!                 IT WILL CONTAIN THE NEWEST ESTIMATE FOR THE SOLUTION
!                 VECTOR.  THIS IS ONLY FOR THE BLACK POINTS.       
!          DB     INPUT D.P. ARRAY.  DB CONTAINS THE VALUE OF THE   
!                 CURRENT PSEUDO-RESIDUAL ON THE BLACK POINTS.      
!          DB1    INPUT/OUTPUT D.P. ARRAY.  DB1 CONTAINS THE PSEUDO-
!                 RESIDUAL ON THE BLACK POINTS FOR THE IN-1 ITERATION 
!                 ON INPUT.  ON OUTPUT, IT IS FOR THE IN+1 ITERATION. 
!          WB     D.P. ARRAY.  WB IS USED FOR COMPUTATIONS INVOLVING
!                 BLACK VECTORS.      
!          TRI    D.P. ARRAY.  STORES THE TRIDIAGONAL MATRIX ASSOCIATED 
!                 WITH CONJUGATE GRADIENT ACCELERATION.   
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(1),JA(1),N,NNB       
      DOUBLE PRECISION A(1),UB(N),UB1(N),DB(NNB),DB1(N),WB(NNB),TRI(2,1)
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER NB,NR,NRP1    
      DOUBLE PRECISION CON,C1,C2,C3,C4,DNRM,GAMOLD,RHOOLD,RHOTMP    
      LOGICAL Q1  
!       
! ... SPECIFICATIONS FOR FUNCTION SUBPROGRAMS   
!       
      DOUBLE PRECISION DDOT 
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN SUBROUTINE RSCG  
!       
! ... COMPUTE NEW ESTIMATE FOR CME IF ADAPT = .TRUE.      
!       
      NB = NNB    
      NR = N-NB   
      NRP1 = NR+1 
      IF (ADAPT) CALL CHGCON (TRI,GAMOLD,RHOOLD,2)
!       
! ... TEST FOR STOPPING     
!       
      DELNNM = DDOT(NB,DB,1,DB,1)     
      DNRM = DELNNM 
      CON = CME   
      CALL PSTOP (NB,UB(NRP1),DNRM,CON,2,Q1)    
      IF (HALT) GO TO 30    
!       
! ... COMPUTE ACCELERATION PARAMETERS 
!       
      CALL VFILL (NR,UB1,0.D0)
      CALL PRSRED (NB,NR,IA,JA,A,DB,UB1)
      CALL VFILL (NB,WB,0.D0) 
      CALL PRSBLK (NB,NR,IA,JA,A,UB1,WB)
      DNRM = DDOT(NB,DB,1,WB,1)       
      IF (ISYM.EQ.0) GO TO 10 
      RHOTMP = DDOT(NB,WB,1,DB1,1)    
      CALL PARCON (DNRM,C1,C2,C3,C4,GAMOLD,RHOTMP,2)      
      RHOOLD = RHOTMP       
      GO TO 20    
   10 CALL PARCON (DNRM,C1,C2,C3,C4,GAMOLD,RHOOLD,2)      
!       
! ... COMPUTE UB(IN+1) AND DB(IN+1)   
!       
   20 CALL SUM3 (NB,C1,DB,C2,UB(NRP1),C3,UB1(NRP1))       
      CALL SUM3 (NB,C1,WB,C4,DB,C3,DB1) 
!       
! ... OUTPUT INTERMEDIATE INFORMATION 
!       
   30 CALL ITERM (NB,A(NRP1),UB(NRP1),WB,6)     
!       
      RETURN      
      END 
      SUBROUTINE ITRSSI (N,NNB,IA,JA,A,RHS,UB,UB1,DB)     
!       
! ... FUNCTION:   
!       
!          THIS SUBROUTINE, ITRSSI, PERFORMS ONE ITERATION OF THE   
!          REDUCED SYSTEM SEMI-ITERATION ALGORITHM.  IT IS
!          CALLED BY RSSI.  
!       
! ... PARAMETER LIST:       
!       
!          N      INPUT INTEGER.  DIMENSION OF THE MATRIX.
!          NB     INPUT INTEGER.  CONTAINS THE NUMBER OF BLACK POINTS 
!                 IN THE RED-BLACK MATRIX. (= NNB)
!          IA,JA  INPUT INTEGER VECTORS.  THE TWO INTEGER ARRAYS OF 
!                 THE SPARSE MATRIX REPRESENTATION.       
!          A      INPUT D.P. VECTOR.  THE D.P. ARRAY OF THE SPARSE  
!                 MATRIX REPRESENTATION.
!          RHS    INPUT D.P. VECTOR.  CONTAINS THE RIGHT HAND SIDE  
!                 OF THE MATRIX PROBLEM.
!          UB     INPUT D.P. VECTOR.  CONTAINS THE ESTIMATE FOR THE 
!                 SOLUTION ON THE BLACK POINTS AFTER IN ITERATIONS. 
!          UB1    INPUT/OUTPUT D.P. VECTOR.  ON INPUT, UB1 CONTAINS THE 
!                 SOLUTION VECTOR AFTER IN-1 ITERATIONS.  ON OUTPUT,
!                 IT WILL CONTAIN THE NEWEST ESTIMATE FOR THE SOLUTION
!                 VECTOR.  THIS IS ONLY FOR THE BLACK POINTS.       
!          DB     INPUT D.P. ARRAY.  DB CONTAINS THE VALUE OF THE   
!                 CURRENT PSEUDO-RESIDUAL ON THE BLACK POINTS.      
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(1),JA(1),N,NNB       
      DOUBLE PRECISION A(1),RHS(N),UB(N),UB1(N),DB(NNB)   
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER NB,NR,NRP1    
      DOUBLE PRECISION CONST,C1,C2,C3,DNRM      
      LOGICAL Q1  
!       
! ... SPECIFICATIONS FOR FUNCTION SUBPROGRAMS   
!       
      DOUBLE PRECISION DDOT 
      LOGICAL TSTCHG
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN SUBROUTINE RSSI  
!       
! ... COMPUTE UR(IN) INTO UB
!       
      NB = NNB    
      NR = N-NB   
      NRP1 = NR+1 
      CALL DCOPY (NR,RHS,1,UB,1)      
      CALL PRSRED (NB,NR,IA,JA,A,UB(NRP1),UB)   
!       
! ... COMPUTE PSEUDO-RESIDUAL, DB(IN) 
!       
      CALL DCOPY (NB,RHS(NRP1),1,DB,1)
      CALL PRSBLK (NB,NR,IA,JA,A,UB,DB) 
      CALL VEVMW (NB,DB,UB(NRP1))     
!       
! ... TEST FOR STOPPING     
!       
      DELNNM = DDOT(NB,DB,1,DB,1)     
      DNRM = DELNNM 
      CONST = CME 
      CALL PSTOP (NB,UB(NRP1),DNRM,CONST,2,Q1)  
      IF (HALT) GO TO 20    
      IF (.NOT.ADAPT) GO TO 10
!       
! ... TEST TO CHANGE PARAMETERS       
!       
      IF (.NOT.TSTCHG(2)) GO TO 10    
!       
! ... CHANGE PARAMETERS     
!       
      CALL VFILL (NR,UB1,0.D0)
      CALL PRSRED (NB,NR,IA,JA,A,DB,UB1)
      DNRM = DDOT(NR,UB1,1,UB1,1)     
      CALL CHGSI (DNRM,2)   
      IF (.NOT.ADAPT) GO TO 10
!       
! ... COMPUTE UB(N+1) AFTER CHANGING PARAMETERS 
!       
      CALL DCOPY (NB,UB(NRP1),1,UB1(NRP1),1)    
      CALL DAXPY (NB,GAMMA,DB,1,UB1(NRP1),1)    
      GO TO 20    
!       
! ... COMPUTE UB(N+1) WITHOUT CHANGE OF PARAMETERS
!       
   10 CALL PARSI (C1,C2,C3,2) 
      CALL SUM3 (NB,C1,DB,C2,UB(NRP1),C3,UB1(NRP1))       
!       
! ... OUTPUT INTERMEDIATE INFORMATION 
!       
   20 CALL ITERM (NB,A(NRP1),UB(NRP1),DB,7)     
!       
      RETURN      
      END 
      INTEGER FUNCTION BISRCH (N,K,L) 
!       
! ... BISRCH IS AN INTEGER FUNCTION WHICH USES A BISECTION SEARCH   
!     TO FIND THE ENTRY J IN THE ARRAY K SUCH THAT THE VALUE L IS   
!     GREATER THAN OR EQUAL TO K(J) AND STRICTLY LESS THAN K(J+1).  
!       
! ... PARAMETER LIST:       
!       
!          N      INTEGER LENGTH OF VECTOR K    
!          K      INTEGER VECTOR      
!          L      INTEGER CONSTANT SUCH THAT  K(J) .GE. L .LT. K(J+1) 
!                 WITH J RETURNED AS VALUE OF INTEGER FUNCTION BISRCH 
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER N,L,K(N)      
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER JLEFT,JMID,JRIGHT       
!       
      JLEFT = 1   
      JRIGHT = N  
      IF (N.EQ.2) GO TO 40  
      JMID = (N+1)/2
!       
   10 IF (L.GE.K(JMID)) GO TO 20      
!       
! ...... L .GE. K(LEFT)  AND  L .LT. K(JMID)    
!       
      JRIGHT = JMID 
      GO TO 30    
!       
! ...... L .GE. K(JMID)  AND  L .LT. K(JRIGHT)  
!       
   20 JLEFT = JMID
!       
! ...... TEST FOR CONVERGENCE 
!       
   30 IF (JRIGHT-JLEFT.EQ.1) GO TO 40 
      JMID = JLEFT+(JRIGHT-JLEFT+1)/2 
      GO TO 10    
!       
! ...... BISECTION SEARCH FINISHED    
!       
   40 BISRCH = JLEFT
!       
      RETURN      
      END 
      DOUBLE PRECISION FUNCTION CHEBY (QA,QT,RRR,IP,CME,SME)
!       
!     COMPUTES THE SOLUTION TO THE CHEBYSHEV EQUATION     
!       
! ... PARAMETER LIST:       
!       
!          QA     RATIO OF PSEUDO-RESIDUALS     
!          QT     VIRTUAL SPECTRAL RADIUS       
!          RRR    ADAPTIVE PARAMETER  
!          IP     NUMBER OF ITERATIONS SINCE LAST CHANGE OF 
!                     PARAMETERS      
!          CME,   ESTIMATES FOR THE LARGEST AND SMALLEST EIGEN-     
!          SME      VALUES OF THE ITERATION MATRIX
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IP  
      DOUBLE PRECISION CME,QA,QT,RRR,SME
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      DOUBLE PRECISION X,Y,Z
!       
      Z = .5D0*(QA+DSQRT(DABS(QA**2-QT**2)))*(1.D0+RRR**IP) 
      X = Z**(1.D0/DBLE(FLOAT(IP)))   
      Y = (X+RRR/X)/(1.D0+RRR)
!       
      CHEBY = .5D0*(CME+SME+Y*(2.D0-CME-SME))   
!       
      RETURN      
      END













      SUBROUTINE CHGCON (TRI,GAMOLD,RHOOLD,IBMTH) 
!       
!     COMPUTES THE NEW ESTIMATE FOR THE LARGEST EIGENVALUE FOR      
!     CONJUGATE GRADIENT ACCELERATION.
!       
! ... PARAMETER LIST:       
!       
!          TRI    TRIDIAGONAL MATRIX ASSOCIATED WITH THE EIGENVALUES
!                    OF THE CONJUGATE GRADIENT POLYNOMIAL 
!          GAMOLD 
!            AND  
!          RHOOLD PREVIOUS VALUES OF ACCELERATION PARAMETERS
!          IBMTH  INDICATOR OF BASIC METHOD BEING ACCELERATED BY CG 
!                      IBMTH = 1,  JACOBI       
!                            = 2,  REDUCED SYSTEM 
!                            = 3,  SSOR 
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IBMTH !1
      DOUBLE PRECISION TRI(2,IN),GAMOLD,RHOOLD
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER IB2,IB3,IER,IP
      DOUBLE PRECISION CMOLD,END,START,EIGVSS,EIGVNS      
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
!       
      GO TO (10,20,30), IBMTH 
!       
! ... JACOBI CONJUGATE GRADIENT       
!       
   10 START = CME !here is for JCG
      IP = IN
!write(6,*) IP,CME,CMOLD,ZETA
      GO TO 40    
!       
! ... REDUCED SYSTEM CG     
!       
   20 START = CME**2
      IP = IN     
      GO TO 40    
!       
! ... SSOR CG     
!       
   30 IF (ADAPT) START = SPR
      IF (.NOT.ADAPT) START = SPECR   
      IP = IN-IS  
!       
! ... DEFINE THE MATRIX     
!       
   40 IF (IP.GE.2) GO TO 60 
      IF (IP.EQ.1) GO TO 50 
!       
! ... IP = 0      
!       
      END = 0.D0  
      CMOLD = 0.D0
      GO TO 110   
!       
! ... IP = 1      
!       
   50 END = 1.D0-1.D0/GAMMA 
      TRI(1,1) = END
      TRI(2,1) = 0.D0       
      GO TO 110   
!       
! ... IP > 1      
!       
   60 IF ((IP.GT.2).AND.(DABS(START-CMOLD).LE.ZETA*START)) GO TO 120
      CMOLD = START
!       
! ... COMPUTE THE LARGEST EIGENVALUE  
!
      TRI(1,IP) = 1.D0-1.D0/GAMMA     
      TRI(2,IP) = (RHO-1.D0)/(RHO*RHOOLD*GAMMA*GAMOLD)    
      IF (ISYM.NE.0) GO TO 80 
      END = EIGVSS(IP,TRI,START,ZETA,ITMAX,IER) 
      IF (IER.EQ.0) GO TO 100 
      IF (LEVEL.GE.2) WRITE (NOUT,70) IER       
   70 FORMAT (/10X,'DIFFICULTY IN COMPUTATION OF MAXIMUM EIGENVALUE'/15X&
     &   ,'OF ITERATION MATRIX'/10X,'SUBROUTINE ZBRENT RETURNED IER =', &
     &   I5)      
      GO TO 100   
   80 IB2 = 1+IP  
      IB3 = IB2+IP/2+1      
      END = EIGVNS(IP,TRI,TRI(1,IB2),TRI(1,IB3),IER)      
      IF (IER.EQ.0) GO TO 100 
      IF (LEVEL.GE.2) WRITE (NOUT,90) IER       
   90 FORMAT (/10X,'DIFFICULTY IN COMPUTATION OF MAXIMUM EIGENVALUE'/15X&
     &   ,'OF ITERATION MATRIX'/10X,'SUBROUTINE EQRT1S RETURNED IER =', &
     &   I5)      
  100 CONTINUE    
      IF (IER.NE.0) GO TO 130 
!       
! ... SET SPECTRAL RADIUS FOR THE VARIOUS METHODS 
!       
  110 IF (IBMTH.EQ.1) CME = END       
      IF (IBMTH.EQ.2) CME = DSQRT(DABS(END))    
      IF (IBMTH.EQ.3.AND.ADAPT) SPR = END       
      IF (IBMTH.EQ.3.AND..NOT.ADAPT) SPECR = END
      RETURN      
!       
! ... RELATIVE CHANGE IN CME IS LESS THAN ZETA.  THEREFORE STOP     
!     CHANGING.   
!       
  120 ADAPT = .FALSE.       
      PARTAD = .FALSE.
      RETURN      
!       
! ... ESTIMATE FOR CME > 1.D0.  THEREFORE NEED TO STOP ADAPTIVE     
!     PROCEDURE AND KEEP OLD VALUE OF CME.      
!       
  130 ADAPT = .FALSE.       
      PARTAD = .FALSE.      
      IF (LEVEL.GE.2) WRITE (NOUT,140) IN,START 
  140 FORMAT (/10X,'ESTIMATE OF MAXIMUM EIGENVALUE OF JACOBI   '/15X, &
     &   'MATRIX (CME) NOT ACCURATE'/10X,       &
     &   'ADAPTIVE PROCEDURE TURNED OFF AT ITERATION ',I5/10X,      &
     &   'FINAL ESTIMATE OF MAXIMUM EIGENVALUE =',D15.7/) 
!       
      RETURN      
      END











      SUBROUTINE CHGSI (DTNRM,IBMTH)  
!       
! ... COMPUTES NEW CHEBYSHEV ACCELERATION PARAMETERS ADAPTIVELY.    
!       
! ... PARAMETER LIST:       
!       
!          DTNRM  NUMERATOR OF RAYLEIGH QUOTIENT
!          IBMTH  INDICATOR OF BASIC METHOD BEING ACCELERATED BY SI 
!                      IBMTH = 1,   JACOBI      
!                            = 2,   REDUCED SYSTEM
!                            = 3,   SYMMETRIC SOR 
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IBMTH 
      DOUBLE PRECISION DTNRM
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      DOUBLE PRECISION CMOLD,ZM1,ZM2  
!       
! ... SPECIFICATIONS FOR FUNCTION SUBPROGRAMS   
!       
      DOUBLE PRECISION CHEBY
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
!       
      GO TO (10,30,50), IBMTH 
!       
!     --------------------- 
! ... JACOBI SEMI-ITERATIVE 
!     --------------------- 
!       
! ... CHEBYSHEV EQUATION    
!       
   10 CONTINUE    
      IF (IN.EQ.0) ZM1 = CME
      IF (IN.NE.0) ZM1 = CHEBY(QA,QT,RRR,IN-IS,CME,SME)   
!       
! ... RAYLEIGH QUOTIENT     
!       
      ZM2 = DTNRM/DELNNM    
!       
! ... COMPUTATION OF ITERATIVE PARAMETERS       
!       
      CMOLD = CME 
      CME = DMAX1(ZM1,ZM2,CMOLD)      
      IF (CME.GE.1.D0) GO TO 20       
      IF (CASEII) SME = -CME
      SIGE = (CME-SME)/(2.D0-CME-SME) 
      GAMMA = 2.D0/(2.D0-CME-SME)     
      RRR = (1.D0-DSQRT(DABS(1.D0-SIGE*SIGE)))/(1.D0+DSQRT(DABS(1.D0- &
     &   SIGE*SIGE)))       
      IS = IN     
      DELSNM = DELNNM       
      RHO = 1.D0  
      IF (LEVEL.GE.2) WRITE (NOUT,90) IN,ZM1,ZM2,CME,GAMMA,CME      
      RETURN      
!       
! ... ADAPTIVE PROCEDURE FAILED FOR JACOBI SI   
!       
   20 CME = CMOLD 
      ADAPT = .FALSE.       
      IF (LEVEL.GE.2) WRITE (NOUT,110) IN,CME   
      RETURN      
!       
!     -----------------------------   
! ... REDUCED SYSTEM SEMI-ITERATIVE   
!     -----------------------------   
!       
! ... CHEBYSHEV EQUATION    
!       
   30 CONTINUE    
      IF (IN.EQ.0) ZM1 = CME
      IF (IN.NE.0) ZM1 = CHEBY(QA,QT,RRR,2*(IN-IS),0.D0,0.D0)       
!       
! ... RAYLEIGH QUOTIENT     
!       
      ZM2 = DSQRT(DABS(DTNRM/DELNNM)) 
!       
! ... COMPUTATION OF NEW ITERATIVE PARAMETERS   
!       
      CMOLD = CME 
      CME = DMAX1(ZM1,ZM2,CMOLD)      
      IF (CME.GE.1.D0) GO TO 40       
      SIGE = CME*CME/(2.D0-CME*CME)   
      GAMMA = 2.D0/(2.D0-CME*CME)     
      RRR = (1.D0-DSQRT(DABS(1.D0-CME*CME)))/(1.D0+DSQRT(DABS(1.D0-CME* &
     &   CME)))   
      IS = IN     
      DELSNM = DELNNM       
      RHO = 1.D0  
      IF (LEVEL.GE.2) WRITE (NOUT,90) IN,ZM1,ZM2,CME,GAMMA,CME      
      RETURN      
!       
! ... ADAPTIVE PROCEDURE FAILED FOR REDUCED SYSTEM SI     
!       
   40 CME = CMOLD 
      ADAPT = .FALSE.       
      IF (LEVEL.GE.2) WRITE (NOUT,110) IN,CME   
      RETURN      
!       
!     -----------------------------   
! ... SYMMETRIC SOR SEMI-ITERATIVE    
!     ----------------------------    
!       
   50 CONTINUE    
      IF (SPECR.EQ.0.D0) SPECR = .171572875D0   
      IF (IN.EQ.0) GO TO 60 
      ZM1 = CHEBY(QA,QT,RRR,IN-IS,SPECR,0.D0)   
      GO TO 70    
   60 ZM1 = SPECR 
      SPR = SPECR 
!       
! ... RAYLEIGH QUOTIENT     
!       
   70 ZM2 = DTNRM/DELNNM    
!       
! ... COMPUTATION OF NEW ESTIMATE FOR SPECTRAL RADIUS     
!       
      IF (ADAPT) GO TO 80   
!       
! ... PARTIALLY ADAPTIVE SSOR SI      
!       
      SPECR = DMAX1(ZM1,ZM2,SPECR)    
      IS = IN+1   
      DELSNM = DELNNM       
      IF (LEVEL.GE.2) WRITE (NOUT,100) IN,ZM1,ZM2,CME,SPECR 
      RETURN      
!       
! ... FULLY ADAPTIVE SSOR SI
!       
   80 SPR = DMAX1(ZM1,ZM2,SPR)
      RETURN      
!       
! ... FORMAT STATEMENTS     
!       
   90 FORMAT (/30X,'PARAMETERS WERE CHANGED AT ITERATION NO.',I5/35X, &
     &   'SOLUTION TO CHEBYSHEV EQN.       =',D15.7/35X,  &
     &   'SOLUTION TO RAYLEIGH QUOTIENT    =',D15.7/35X,  &
     &   'NEW ESTIMATE FOR CME             =',D15.7/35X,  &
     &   'NEW ESTIMATE FOR GAMMA           =',D15.7/35X,  &
     &   'NEW ESTIMATE FOR SPECTRAL RADIUS =',D15.7/)     
!       
  100 FORMAT (/30X,'PARAMETERS WERE CHANGED AT ITERATION NO.',I5/35X, &
     &   'SOLUTION TO CHEBYSHEV EQN.       =',D15.7/35X,  &
     &   'SOLUTION TO RAYLEIGH QUOTIENT    =',D15.7/35X,  &
     &   'NEW ESTIMATE FOR CME             =',D15.7/35X,  &
     &   'NEW ESTIMATE FOR SPECTRAL RADIUS =',D15.7/)     
!       
  110 FORMAT (/10X,'ESTIMATE OF MAXIMUM EIGENVALUE OF JACOBI   '/15X, &
     &   'MATRIX (CME) TOO LARGE'/10X,&
     &   'ADAPTIVE PROCEDURE TURNED OFF AT ITERATION ',I5/10X,      &
     &   'FINAL ESTIMATE OF MAXIMUM EIGENVALUE =',D15.7/) 
!       
      END 
      LOGICAL FUNCTION CHGSME (OLDNRM,ICNT)     
!       
! ... THIS FUNCTION TESTS FOR JACOBI SI WHETHER SME SHOULD BE CHANGED 
! ... WHEN CASEII = .FALSE..  IF THE TEST IS POSITIVE THE NEW VALUE 
! ... OF SME IS COMPUTED.   
!       
! ... PARAMETER LIST:       
!       
!          OLDNRM SQUARE OF THE NORM OF THE PSEUDO-RESIDUAL 
!                    AT THE LAST ITERATION      
!          ICNT   NUMBER OF ITERATIONS SINCE LAST CHANGE OF 
!                    PARAMETERS       
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER ICNT
      DOUBLE PRECISION OLDNRM 
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER IP  
      DOUBLE PRECISION Q,RN,SM1,SM2,WP,Z
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
!       
      CHGSME = .FALSE.      
      RN = DSQRT(DELNNM/OLDNRM)       
      IF (.NOT.(QA.GT.1.D0.AND.RN.GT.1.D0)) RETURN
      IF (IN.LE.IS+2) RETURN
!       
      ICNT = ICNT+1 
      IF (ICNT.LT.3) RETURN 
!       
! ... CHANGE SME IN J-SI ADAPTIVE PROCEDURE     
!       
      CHGSME = .TRUE.       
      SM1 = 0.D0  
      SM2 = 0.D0  
      IF (SME.GE.CME) GO TO 10
!       
! ... COMPUTE SM1 
!       
      IP = IN-IS  
      Q = QA*(1.D0+RRR**IP)/(2.D0*DSQRT(RRR**IP)) 
      Z = (Q+DSQRT(Q**2-1.D0))**(1.D0/DBLE(FLOAT(IP)))    
      WP = (Z**2+1.D0)/(2.D0*Z)       
      SM1 = .5D0*(CME+SME-WP*(CME-SME)) 
!       
! ... COMPUTE SM2 
!       
      Q = RN*(1.D0+RRR**IP)/((1.D0+RRR**(IP-1))*DSQRT(RRR)) 
      WP = (Q**2+1.D0)/(2.D0*Q)       
      SM2 = .5D0*(CME+SME-WP*(CME-SME)) 
!       
   10 SME = DMIN1(1.25D0*SM1,1.25D0*SM2,SME,-1.D0)
      SIGE = (CME-SME)/(2.D0-CME-SME) 
      GAMMA = 2.D0/(2.D0-CME-SME)     
      RRR = (1.D0-DSQRT(1.D0-SIGE**2))/(1.D0+DSQRT(1.D0-SIGE**2))   
      IS = IN     
      DELSNM = DELNNM       
      RHO = 1.D0  
!       
      IF (LEVEL.GE.2) WRITE (NOUT,20) IN,SM1,SM2,SME      
!       
   20 FORMAT (/30X,'ESTIMATE OF SMALLEST EIGENVALUE OF JACOBI'/37X, &
     &   'MATRIX (SME) CHANGED AT ITERATION ',I5/35X,     &
     &   'FIRST ESTIMATE OF SME            =',D15.7/35X,  &
     &   'SECOND ESTIMATE OF SME           =',D15.7/35X,  &
     &   'NEW ESTIMATE OF SME              =',D15.7/)     
!       
      RETURN      
      END









      SUBROUTINE DAXPY (N,DA,DX,INCX,DY,INCY)   
!       
!     OVERWRITE DOUBLE PRECISION DY WITH DOUBLE PRECISION DA*DX + DY. 
!       
      DOUBLE PRECISION DX(1),DY(1),DA 
      IF (N.LE.0.OR.DA.EQ.0.D0) RETURN
      IF (INCX.EQ.INCY) IF (INCX-1) 10 , 30 , 70
   10 CONTINUE    
!       
!        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.     
!       
      IX = 1      
      IY = 1      
      IF (INCX.LT.0) IX = (-N+1)*INCX+1 
      IF (INCY.LT.0) IY = (-N+1)*INCY+1 
      DO 20 I = 1,N 
         DY(IY) = DY(IY)+DA*DX(IX)    
         IX = IX+INCX       
         IY = IY+INCY       
   20 CONTINUE    
      RETURN      
!       
!        CODE FOR BOTH INCREMENTS EQUAL TO 1    
!       
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4. 
!       
   30 M = N-(N/4)*4 
      IF (M.EQ.0) GO TO 50  
      DO 40 I = 1,M 
         DY(I) = DY(I)+DA*DX(I)       
   40 CONTINUE    
      IF (N.LT.4) RETURN    
   50 MP1 = M+1   
      DO 60 I = MP1,N,4     
         DY(I) = DY(I)+DA*DX(I)       
         DY(I+1) = DY(I+1)+DA*DX(I+1) 
         DY(I+2) = DY(I+2)+DA*DX(I+2) 
         DY(I+3) = DY(I+3)+DA*DX(I+3) 
   60 CONTINUE    
      RETURN      
!       
!        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.    
!       
   70 CONTINUE    
      NS = N*INCX 
      DO 80 I = 1,NS,INCX   
         DY(I) = DA*DX(I)+DY(I)       
   80 CONTINUE    
      RETURN      
      END









      SUBROUTINE DCOPY (N,DX,INCX,DY,INCY)      
!       
!     COPY DOUBLE PRECISION DX TO DOUBLE PRECISION DY.    
!       
      DOUBLE PRECISION DX(1+N),DY(1+N)
!     DOUBLE PRECISION DX(1),DY(1)  !old version

      IF (N.LE.0) RETURN    
      IF (INCX.EQ.INCY) IF (INCX-1) 10 , 30 , 70
   10 CONTINUE    
!       
!        CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.      
!       
      IX = 1      
      IY = 1      
      IF (INCX.LT.0) IX = (-N+1)*INCX+1 
      IF (INCY.LT.0) IY = (-N+1)*INCY+1 
      DO 20 I = 1,N 
         DY(IY) = DX(IX)    
         IX = IX+INCX       
         IY = IY+INCY       
   20 CONTINUE    
      RETURN      
!       
!        CODE FOR BOTH INCREMENTS EQUAL TO 1    
!       
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7. 
!       
   30 M = N-(N/7)*7 
      IF (M.EQ.0) GO TO 50  
      DO 40 I = 1,M 
         DY(I) = DX(I)      
   40 CONTINUE    
      IF (N.LT.7) RETURN    
   50 MP1 = M+1   
      DO 60 I = MP1,N,7     
         DY(I) = DX(I)      
         DY(I+1) = DX(I+1)  
         DY(I+2) = DX(I+2)  
         DY(I+3) = DX(I+3)  
         DY(I+4) = DX(I+4)  
         DY(I+5) = DX(I+5)  
         DY(I+6) = DX(I+6)  
   60 CONTINUE    
      RETURN      
!       
!        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.    
!       
   70 CONTINUE    
      NS = N*INCX 
      DO 80 I = 1,NS,INCX   
         DY(I) = DX(I)      
   80 CONTINUE    
      RETURN      
      END 
      DOUBLE PRECISION FUNCTION DDOT (N,DX,INCX,DY,INCY)  
!       
!     RETURNS THE DOT PRODUCT OF DOUBLE PRECISION DX AND DY.
!       
      DOUBLE PRECISION DX(N+5),DY(N+5)
      DDOT = 0.D0 
      IF (N.LE.0) RETURN    
      IF (INCX.EQ.INCY) IF (INCX-1) 10 , 30 , 70
   10 CONTINUE    
!       
!         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.     
!       
      IX = 1      
      IY = 1      
      IF (INCX.LT.0) IX = (-N+1)*INCX+1 
      IF (INCY.LT.0) IY = (-N+1)*INCY+1 
      DO 20 I = 1,N 
         DDOT = DDOT+DX(IX)*DY(IY)    
         IX = IX+INCX       
         IY = IY+INCY       
   20 CONTINUE    
      RETURN      
!       
!        CODE FOR BOTH INCREMENTS EQUAL TO 1.   
!       
!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5. 
!       
   30 M = N-(N/5)*5 
      IF (M.EQ.0) GO TO 50  
      DO 40 I = 1,M 
         DDOT = DDOT+DX(I)*DY(I)      
   40 CONTINUE    
      IF (N.LT.5) RETURN    
   50 MP1 = M+1   
      DO 60 I = MP1,N,5     
         DDOT = DDOT+DX(I)*DY(I)+DX(I+1)*DY(I+1)+DX(I+2)*DY(I+2)+DX(I+3)&
     &      *DY(I+3)+DX(I+4)*DY(I+4)  
   60 CONTINUE    
      RETURN      
!       
!         CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.       
!       
   70 CONTINUE    
      NS = N*INCX 
      DO 80 I = 1,NS,INCX   
         DDOT = DDOT+DX(I)*DY(I)      
   80 CONTINUE    
      RETURN      
      END 
      DOUBLE PRECISION FUNCTION DETERM (N,TRI,XLMDA)      
!       
!     THIS SUBROUTINE COMPUTES THE DETERMINANT OF A SYMMETRIC       
!     TRIDIAGONAL MATRIX GIVEN BY TRI. DET(TRI - XLMDA*I) = 0       
!       
! ... PARAMETER LIST
!       
!          N      ORDER OF TRIDIAGONAL SYSTEM   
!          TRI    SYMMETRIC TRIDIAGONAL MATRIX OF ORDER N 
!          XLMDA  ARGUMENT FOR CHARACTERISTIC EQUATION    
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER N   
      DOUBLE PRECISION TRI(2,N),XLMDA
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER ICNT,L,NM1    
      DOUBLE PRECISION D1,D2,D3       
!       
      NM1 = N-1   
      D2 = TRI(1,N)-XLMDA   
      D1 = D2*(TRI(1,NM1)-XLMDA)-TRI(2,N)       
      IF (N.EQ.2) GO TO 20  
!       
! ... BEGINNING OF LOOP     
!       
      DO 10 ICNT = 2,NM1    
         L = NM1-ICNT+2     
         D3 = D2  
         D2 = D1  
         D1 = (TRI(1,L-1)-XLMDA)*D2-D3*TRI(2,L) 
   10 CONTINUE    
!       
! ... DETERMINANT COMPUTED  
!       
   20 DETERM = D1 
!       
      RETURN      
      END













      SUBROUTINE DFAULT (IPARM,RPARM) 
!       
! ... THIS SUBROUTINE SETS THE DEFAULT VALUES OF IPARM AND RPARM.   
!       
! ... PARAMETER LIST:       
!       
!          IPARM  
!           AND   
!          RPARM  ARRAYS SPECIFYING OPTIONS AND TOLERANCES
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IPARM(12)     
      DOUBLE PRECISION RPARM(12)      
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
!       
!     DRELPR  - COMPUTER PRECISION (APPROX.)    
!     IF INSTALLER OF PACKAGE DOES NOT KNOW DRELPR VALUE, 
!     AN APPROXIMATE VALUE CAN BE DETERMINED FROM A SIMPLE
!     FORTRAN PROGRAM SUCH AS 
!       
!     DOUBLE PRECISION DRELPR, TEMP   
!     DRELPR = 1.0D0
!   2 DRELPR = 0.5D0*DRELPR 
!     TEMP = DRELPR + 1.0D0 
!     IF(TEMP .GT. 1.0D0)  GO TO 2    
!     WRITE(6,3) DRELPR     
!   3 FORMAT(5X,D15.8)      
!     STOP
!     END 
!       
!     SOME VALUES ARE:      
!       
!     DRELPR = 1.26D-29  FOR CDC CYBER 170/750  (APPROX.) 2**-96    
!            = 2.22D-16  FOR DEC 10             (APPROX.) 2**-52    
!            = 7.11D-15  FOR VAX 11/780         (APPROX.) 2**-47    
!            = 1.14D-13  FOR IBM 370/158        (APPROX.) 2**-43    
!       
!             *** SHOULD BE CHANGED FOR OTHER MACHINES ***
!       
!     TO FACILITATE CONVERGENCE, RPARM(1) SHOULD BE SET TO
!          500.*DRELPR OR LARGER      
!       
      DRELPR = 0.11102230D-15   !7.11D-15
!       
      IPARM(1) = 100
      IPARM(2) = 0
      IPARM(3) = 0
      IPARM(4) = 6
      IPARM(5) = 0
      IPARM(6) = 1
      IPARM(7) = 1
      IPARM(8) = 0
      IPARM(9) = -1 
      IPARM(10) = 0 
      IPARM(11) = 0 
      IPARM(12) = 0 
!       
      RPARM(1) = 0.5D-5     
      RPARM(2) = 0.D0       
      RPARM(3) = 0.D0       
      RPARM(4) = .75D0      
      RPARM(5) = 1.D0       
      RPARM(6) = 0.D0       
      RPARM(7) = .25D0      
      RPARM(8) = 1.D2*DRELPR
      RPARM(9) = 0.D0       
      RPARM(10) = 0.D0      
      RPARM(11) = 0.D0      
      RPARM(12) = 0.D0      
!       
      RETURN      
      END










      SUBROUTINE ECHALL (NN,IA,JA,A,RHS,IPARM,RPARM,ICALL)
!       
! ... THIS ROUTINE INITIALIZES THE ITPACK COMMON BLOCKS FROM THE    
! ... INFORMATION CONTAINED IN IPARM AND RPARM. ECHALL ALSO PRINTS THE
! ... VALUES OF ALL THE PARAMETERS IN IPARM AND RPARM.    
!       
! ... PARAMETER LIST:       
!       
!          IPARM  
!           AND   
!          RPARM  ARRAYS OF PARAMETERS SPECIFYING OPTIONS AND       
!                    TOLERANCES       
!          ICALL  INDICATOR OF WHICH PARAMETERS ARE BEING PRINTED   
!                    ICALL = 1,  INITIAL PARAMETERS       
!                    ICALL = 2,  FINAL PARAMETERS 
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(1),JA(1),IPARM(12),NN,ICALL    
      DOUBLE PRECISION A(1),RHS(NN),RPARM(12)   
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER I,N,NP1,NZRO  
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
!       
      IF (ICALL.NE.1) GO TO 100       
      N = NN      
      NP1 = N+1   
      NZRO = IA(NP1)-1      
!       
! ... INITIALIZE ITPACK COMMON
!       
      ZETA = RPARM(1)       
      CME = RPARM(2)
      SME = RPARM(3)
      FF = RPARM(4) 
      OMEGA = RPARM(5)      
      SPECR = RPARM(6)      
      BETAB = RPARM(7)      
      ITMAX = IPARM(1)      
      LEVEL = IPARM(2)      
      ISYM = IPARM(5)       
!       
      ADAPT = .FALSE.       
      PARTAD = .FALSE.      
      BETADT = .FALSE.      
      IF (IPARM(6).EQ.1.OR.IPARM(6).EQ.3) ADAPT = .TRUE.  
      IF (IPARM(6).EQ.1) BETADT = .TRUE.
      IF (IPARM(6).EQ.2) PARTAD = .TRUE.
!       
      CASEII = .FALSE.      
      IF (IPARM(7).EQ.2) CASEII = .TRUE.
      IF (CASEII) SME = -CME
      IF (.NOT.CASEII.AND.SME.EQ.0.D0) SME = -1.D0
      SPR = SME   
!       
! ... SET REST OF COMMON VARIABLES TO ZERO      
!       
      IN = 0      
      IS = 0      
      HALT = .FALSE.
      BDELNM = 0.D0 
      DELNNM = 0.D0 
      DELSNM = 0.D0 
      GAMMA = 0.D0
      QA = 0.D0   
      QT = 0.D0   
      RHO = 0.D0  
      RRR = 0.D0  
      SIGE = 0.D0 
      STPTST = 0.D0 
      UDNM = 0.D0 
!       
      IF (LEVEL.LE.4) GO TO 80
!       
!     THIS SECTION OF ECHALL CAUSES PRINTING OF THE LINEAR SYSTEM AND 
!     THE ITERATIVE PARAMETERS
!       
      WRITE (NOUT,10)       
   10 FORMAT (///30X,'THE LINEAR SYSTEM IS AS FOLLOWS')   
      WRITE (NOUT,20)       
   20 FORMAT (/2X,'IA ARRAY') 
      WRITE (NOUT,30) (IA(I),I=1,NP1) 
   30 FORMAT (2X,10(2X,I8)) 
      WRITE (NOUT,40)       
   40 FORMAT (/2X,'JA ARRAY') 
      WRITE (NOUT,30) (JA(I),I=1,NZRO)
      WRITE (NOUT,50)       
   50 FORMAT (/2X,' A ARRAY') 
      WRITE (NOUT,60) (A(I),I=1,NZRO) 
   60 FORMAT (2X,5(2X,D20.13))
      WRITE (NOUT,70)       
   70 FORMAT (/2X,'RHS ARRAY')
      WRITE (NOUT,60) (RHS(I),I=1,N)  
   80 WRITE (NOUT,90)       
   90 FORMAT (///30X,'INITIAL ITERATIVE PARAMETERS')      
      GO TO 120   
  100 WRITE (NOUT,110)      
  110 FORMAT (///30X,'FINAL ITERATIVE PARAMETERS')
  120 WRITE (NOUT,130) IPARM(1),LEVEL,IPARM(3),NOUT,ISYM,IPARM(6)   
  130 FORMAT (35X,'IPARM(1)  =',I15,4X,'(ITMAX)'/35X,'IPARM(2)  =',I15, &
     &   4X,'(LEVEL) '/35X,'IPARM(3)  =',I15,4X,'(IRESET)'/35X,     &
     &   'IPARM(4)  =',I15,4X,'(NOUT)  '/35X,'IPARM(5)  =',I15,4X,  &
     &   '(ISYM)  '/35X,'IPARM(6)  =',I15,4X,'(IADAPT)')  
      WRITE (NOUT,140) IPARM(7),IPARM(8),IPARM(9),IPARM(10),IPARM(11),&
     &   IPARM(12)
  140 FORMAT (35X,'IPARM(7)  =',I15,4X,'(ICASE)'/35X,'IPARM(8)  =',I15, &
     &   4X,'(NWKSP)'/35X,'IPARM(9)  =',I15,4X,'(NB)    '/35X,      &
     &   'IPARM(10) =',I15,4X,'(IREMOVE)'/35X,'IPARM(11) =',I15,4X, &
     &   '(ITIME)'/35X,'IPARM(12) =',I15,4X,'(IDGTS)')    
      WRITE (NOUT,150) ZETA,CME,SME,FF,OMEGA,SPECR
  150 FORMAT (35X,'RPARM(1)  =',D15.8,4X,'(ZETA)  '/35X,'RPARM(2)  =',&
     &   D15.8,4X,'(CME)   '/35X,'RPARM(3)  =',D15.8,4X,'(SME)   '/35X, &
     &   'RPARM(4)  =',D15.8,4X,'(FF)    '/35X,'RPARM(5)  =',D15.8,4X,&
     &   '(OMEGA) '/35X,'RPARM(6)  =',D15.8,4X,'(SPECR) ')
      WRITE (NOUT,160) BETAB,RPARM(8),RPARM(9),RPARM(10),RPARM(11), &
     &   RPARM(12)
  160 FORMAT (35X,'RPARM(7)  =',D15.8,4X,'(BETAB) '/35X,'RPARM(8)  =',&
     &   D15.8,4X,'(TOL)'/35X,'RPARM(9)  =',D15.8,4X,'(TIME1)'/35X, &
     &   'RPARM(10) =',D15.8,4X,'(TIME2)'/35X,'RPARM(11) =',D15.8,4X, &
     &   '(DIGIT1)'/35X,'RPARM(12) =',D15.8,4X,'(DIGIT2)')
!       
      RETURN      
      END









      SUBROUTINE ECHOUT (IPARM,RPARM,IMTHD)     
!       
!     THIS ROUTINE INITIALIZES THE ITPACK COMMON BLOCKS FROM THE    
!     INFORMATION CONTAINED IN IPARM AND RPARM. 
!       
! ... PARAMETER LIST:       
!       
!          IPARM  
!           AND   
!          RPARM  ARRAYS OF PARAMETERS SPECIFYING OPTIONS AND       
!                    TOLERANCES       
!          IMTHD  INDICATOR OF METHOD 
!                    IMTHD = 1,  JCG  
!                    IMTHD = 2,  JSI  
!                    IMTHD = 3,  SOR  
!                    IMTHD = 4,  SSORCG 
!                    IMTHD = 5,  SSORSI 
!                    IMTHD = 6,  RSCG 
!                    IMTHD = 7,  RSSI 
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IPARM(12),IMTHD 
      DOUBLE PRECISION RPARM(12)      
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
!       
! ... INITIALIZE ITPACK COMMON
!       
      ZETA = RPARM(1)       
      CME = RPARM(2)
      SME = RPARM(3)
      FF = RPARM(4) 
      OMEGA = RPARM(5)      
      SPECR = RPARM(6)      
      BETAB = RPARM(7)      
      ITMAX = IPARM(1)      
      LEVEL = IPARM(2)      
      ISYM = IPARM(5)       
!       
      ADAPT = .FALSE.       
      PARTAD = .FALSE.      
      BETADT = .FALSE.      
      IF (IPARM(6).EQ.1.OR.IPARM(6).EQ.3) ADAPT = .TRUE.  
      IF (IPARM(6).EQ.1) BETADT = .TRUE.
      IF (IPARM(6).EQ.2) PARTAD = .TRUE.
!       
      CASEII = .FALSE.      
      IF (IPARM(7).EQ.2) CASEII = .TRUE.
      IF (CASEII) SME = -CME
      IF (.NOT.CASEII.AND.SME.EQ.0.D0) SME = -1.D0
      SPR = SME   
!       
! ... SET REST OF COMMON VARIABLES TO ZERO      
!       
      IN = 0      
      IS = 0      
      HALT = .FALSE.
      BDELNM = 0.D0 
      DELNNM = 0.D0 
      DELSNM = 0.D0 
      GAMMA = 0.D0
      QA = 0.D0   
      QT = 0.D0   
      RHO = 0.D0  
      RRR = 0.D0  
      SIGE = 0.D0 
      STPTST = 0.D0 
      UDNM = 0.D0 
      IF (LEVEL.LE.2) RETURN
!       
! ... THIS SECTION OF ECHOUT ECHOES THE INPUT VALUES FOR THE INITIAL
!     ITERATIVE PARAMETERS  
!       
      WRITE (NOUT,10) ISYM,ITMAX,ZETA,ADAPT,CASEII
   10 FORMAT (///30X,'INITIAL ITERATIVE PARAMETERS',3X,   &
     &   'RELEVANT SWITCHES'/35X,'ISYM   =',I15,8X,'IPARM(5)'/35X,  &
     &   'ITMAX  =',I15,8X,'IPARM(1)'/35X,'ZETA   =',D15.8,8X,'RPARM(1)'&
     &   /35X,'ADAPT  =',L15,8X,'IPARM(6)'/35X,'CASEII =',L15,8X,   &
     &   'IPARM(7)')
      GO TO (80,20,100,60,40,80,20), IMTHD      
!       
! ... JSI, RSSI   
!       
   20 WRITE (NOUT,30) FF,CME,SME      
   30 FORMAT (35X,'FF     =',D15.8,8X,'RPARM(4)'/35X,'CME    =',D15.8,8X&
     &   ,'RPARM(2)'/35X,'SME    =',D15.8,8X,'RPARM(3)'///) 
      RETURN      
!       
! ... SSORSI      
!       
   40 WRITE (NOUT,50) PARTAD,FF,CME,OMEGA,SPECR,BETAB,BETADT
   50 FORMAT (35X,'PARTAD =',L15,8X,'IPARM(6)'/35X,'FF     =',D15.8,8X, &
     &   'RPARM(4)'/35X,'CME    =',D15.8,8X,'RPARM(2)'/35X,'OMEGA  =',&
     &   D15.8,8X,'RPARM(5)'/35X,'SPECR  =',D15.8,8X,'RPARM(6)'/35X,&
     &   'BETAB  =',D15.8,8X,'RPARM(7)'/35X,'BETADT =',L15,8X,'IPARM(6)'&
     &   ///)     
      RETURN      
!       
! ... SSORCG      
!       
   60 WRITE (NOUT,70) PARTAD,CME,OMEGA,SPECR,BETAB,BETADT 
   70 FORMAT (35X,'PARTAD =',L15,8X,'IPARM(6)'/35X,'CME    =',D15.8,8X, &
     &   'RPARM(2)'/35X,'OMEGA  =',D15.8,8X,'RPARM(5)'/35X,'SPECR  =',&
     &   D15.8,8X,'RPARM(6)'/35X,'BETAB  =',D15.8,8X,'RPARM(7)'/35X,&
     &   'BETADT =',L15,8X,'IPARM(6)'///)       
      RETURN      
!       
! ... JCG, RSCG   
!       
   80 IF (ADAPT) RETURN     
      WRITE (NOUT,90) CME   
   90 FORMAT (35X,'CME    =',D15.8,8X,'RPARM(2)'///)      
!       
  100 CONTINUE    
      RETURN      
      END 
      DOUBLE PRECISION FUNCTION EIGVNS (N,TRI,D,E2,IER)   
!       
!     COMPUTES THE LARGEST EIGENVALUE OF A SYMMETRIC TRIDIAGONAL MATRIX 
!     FOR CONJUGATE GRADIENT ACCELERATION.      
!       
! ... PARAMETER LIST:       
!       
!          N      ORDER OF TRIDIAGONAL SYSTEM   
!          TRI    SYMMETRIC TRIDIAGONAL MATRIX OF ORDER N 
!          D      ARRAY FOR EQRT1S (NEGATIVE DIAGONAL ELEMENTS)     
!          E2     ARRAY FOR EQRT1S (SUPER DIAGONAL ELEMENTS)
!          IER    ERROR FLAG: ON RETURN, IER=0 INDICATES THAT       
!                    THE LARGEST EIGENVALUE OF TRI WAS FOUND.       
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER N,IER 
      DOUBLE PRECISION TRI(2,1),D(N),E2(N)      
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER I   
!       
      EIGVNS = 0.D0 
!       
      D(1) = -TRI(1,1)      
      DO 10 I = 2,N 
         D(I) = -TRI(1,I)   
         E2(I) = DABS(TRI(2,I))       
   10 CONTINUE    
!       
      CALL EQRT1S (D,E2,N,1,0,IER)    
      EIGVNS = -D(1)
!       
      RETURN      
      END 
      DOUBLE PRECISION FUNCTION EIGVSS (N,TRI,START,ZETA,ITMAX,IER) 
!       
!     COMPUTES THE LARGEST EIGENVALUE OF A SYMMETRIC TRIDIAGONAL MATRIX 
!     FOR CONJUGATE GRADIENT ACCELERATION.      
!     MODIFIED IMSL ROUTINE ZBRENT USED.
!       
! ... PARAMETER LIST:       
!       
!          N      ORDER OF TRIDIAGONAL SYSTEM   
!          TRI    SYMMETRIC TRIDIAGONAL MATRIX OF ORDER N 
!          START  INITIAL LOWER BOUND OF INTERVAL CONTAINING ROOT   
!          ZETA   STOPPING CRITERIA   
!          IER    ERROR FLAG: ON RETURN, IER=0 INDICATES THAT       
!                    THE LARGEST EIGENVALUE OF TRI WAS FOUND.       
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER N,ITMAX,IER   
      DOUBLE PRECISION TRI(2,1),START,ZETA,A,B,EPS
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER MAXFN,NSIG,ITMP 
!       
      EIGVSS = 0.D0 
      ITMP = IFIX(SNGL(-DLOG10(DABS(ZETA))))    
      NSIG = MAX0(ITMP,4)   
      MAXFN = MAX0(ITMAX,50)
!       
!     EPS = DMIN1(ZETA,0.5D-4)
!       
      EPS = 0.0D0 
      A = START   
      B = 1.0D0   
      CALL ZBRENT (N,TRI,EPS,NSIG,A,B,MAXFN,IER)
      EIGVSS = B  
!       
      RETURN      
      END









      SUBROUTINE EQRT1S (D,E2,NN,M,ISW,IERR)    
!       
!   MODIFIED IMSL ROUTINE NAME   - EQRT1S       
!       
!-----------------------------------------------------------------------
!       
!   COMPUTER            - CDC/SINGLE  
!       
!   LATEST REVISION     - JUNE 1, 1980
!       
!   PURPOSE             - SMALLEST OR LARGEST M EIGENVALUES OF A    
!                           SYMMETRIC TRIDIAGONAL MATRIX  
!       
!   USAGE               - CALL EQRT1S (D,E2,N,M,ISW,IER)  
!       
!   ARGUMENTS    D      - INPUT VECTOR OF LENGTH N CONTAINING       
!                           THE DIAGONAL ELEMENTS OF THE MATRIX.  THE 
!                           COMPUTED EIGENVALUES REPLACE THE FIRST M
!                           COMPONENTS OF THE VECTOR D IN NON-      
!                           DECREASING SEQUENCE, WHILE THE REMAINING
!                           COMPONENTS ARE LOST.
!                E2     - INPUT VECTOR OF LENGTH N CONTAINING       
!                           THE SQUARES OF THE OFF-DIAGONAL ELEMENTS
!                           OF THE MATRIX.  INPUT E2 IS DESTROYED.  
!                N      - INPUT SCALAR CONTAINING THE ORDER OF THE  
!                           MATRIX. (= NN)      
!                M      - INPUT SCALAR CONTAINING THE NUMBER OF     
!                           SMALLEST EIGENVALUES DESIRED (M IS      
!                           LESS THAN OR EQUAL TO N).     
!                ISW    - INPUT SCALAR MEANING AS FOLLOWS - 
!                           ISW=1 MEANS THAT THE MATRIX IS KNOWN TO BE
!                             POSITIVE DEFINITE.
!                           ISW=0 MEANS THAT THE MATRIX IS NOT KNOWN
!                             TO BE POSITIVE DEFINITE.    
!                IER    - ERROR PARAMETER. (OUTPUT) (= IERR)
!                           WARNING ERROR       
!                             IER = 601 INDICATES THAT SUCCESSIVE   
!                               ITERATES TO THE K-TH EIGENVALUE WERE NOT
!                               MONOTONE INCREASING. THE VALUE K IS 
!                               STORED IN E2(1).
!                           TERMINAL ERROR      
!                             IER = 602 INDICATES THAT ISW=1 BUT MATRIX 
!                               IS NOT POSITIVE DEFINITE  
!       
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32 
!                       - SINGLE/H36,H48,H60    
!       
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND       
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL  
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!       
!   REMARKS      AS WRITTEN, THE ROUTINE COMPUTES THE M SMALLEST    
!                EIGENVALUES. TO COMPUTE THE M LARGEST EIGENVALUES, 
!                REVERSE THE SIGN OF EACH ELEMENT OF D BEFORE AND   
!                AFTER CALLING THE ROUTINE. IN THIS CASE, ISW MUST  
!                EQUAL ZERO.
!       
!   COPYRIGHT           - 1980 BY IMSL, INC. ALL RIGHTS RESERVED.   
!       
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.    
!       
!-----------------------------------------------------------------------
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
!       
!                                  SPECIFICATIONS FOR ARGUMENTS     
!       
      INTEGER NN,M,ISW,IERR 
      DOUBLE PRECISION D(NN),E2(NN)   
!       
!                                  SPECIFICATIONS FOR LOCAL VARIABLES 
!       
      INTEGER II,I,JJ,J,K1,K,N,IER    
      DOUBLE PRECISION DELTA,DLAM,EP,ERR,F,P,QP,Q,R,S,TOT 
!       
!                                  DRELPR = MACHINE PRECISION       
!                                  FIRST EXECUTABLE STATEMENT       
!       
      N = NN      
      IER = 0     
      DLAM = 0.0D0
      ERR = 0.0D0 
      S = 0.0D0   
!       
!                                  LOOK FOR SMALL SUB-DIAGONAL ENTRIES
!                                  DEFINE INITIAL SHIFT FROM LOWER  
!                                  GERSCHGORIN BOUND.     
!       
      TOT = D(1)  
      Q = 0.0D0   
      J = 0       
      DO 30 I = 1,N 
         P = Q    
         IF (I.EQ.1) GO TO 10 
         IF (P.GT.DRELPR*(DABS(D(I))+DABS(D(I-1)))) GO TO 20
   10    E2(I) = 0.0D0      
!       
!                                  COUNT IF E2(I) HAS UNDERFLOWED   
!       
   20    IF (E2(I).EQ.0.D0) J = J+1   
         Q = 0.0D0
         IF (I.NE.N) Q = DSQRT(DABS(E2(I+1)))   
         TOT = DMIN1(D(I)-P-Q,TOT)    
   30 CONTINUE    
      IF (ISW.EQ.1.AND.TOT.LT.0.0D0) GO TO 50   
      DO 40 I = 1,N 
         D(I) = D(I)-TOT    
   40 CONTINUE    
      GO TO 60    
   50 TOT = 0.0D0 
   60 DO 200 K = 1,M
!       
!                                  NEXT QR TRANSFORMATION 
!       
   70    TOT = TOT+S
         DELTA = D(N)-S     
         I = N    
         F = DABS(DRELPR*TOT) 
         IF (DLAM.LT.F) DLAM = F      
         IF (DELTA.GT.DLAM) GO TO 90  
         IF (DELTA.GE.(-DLAM)) GO TO 170
         IER = 602
         IF (LEVEL.GE.1) WRITE (NOUT,80)
   80    FORMAT ('0','*** W A R N I N G ************'/' ',&
     &      '    IN ITPACK ROUTINE EQRT1S  '/' ', &
     &      '    PARAMETER ISW = 1 BUT MATRIX   '/' ',    &
     &      '    NOT POSITIVE DEFINITE')
         GO TO 210
!       
!                                  REPLACE SMALL SUB-DIAGONAL SQUARES 
!                                  BY ZERO TO REDUCE THE INCIDENCE OF 
!                                  UNDERFLOWS   
!       
   90    IF (K.EQ.N) GO TO 110
         K1 = K+1 
         DO 100 J = K1,N    
            IF (E2(J).LE.(DRELPR*(D(J)+D(J-1)))**2) E2(J) = 0.0D0   
  100    CONTINUE 
  110    F = E2(N)/DELTA    
         QP = DELTA+F       
         P = 1.0D0
         IF (K.EQ.N) GO TO 140
         K1 = N-K 
         DO 130 II = 1,K1   
            I = N-II
            Q = D(I)-S-F    
            R = Q/QP
            P = P*R+1.0D0   
            EP = F*R
            D(I+1) = QP+EP  
            DELTA = Q-EP    
            IF (DELTA.GT.DLAM) GO TO 120
            IF (DELTA.GE.(-DLAM)) GO TO 170     
            IER = 602       
            IF (LEVEL.GE.0) WRITE (NOUT,80)     
            GO TO 210       
  120       F = E2(I)/Q     
            QP = DELTA+F    
            E2(I+1) = QP*EP 
  130    CONTINUE 
  140    D(K) = QP
         S = QP/P 
         IF (TOT+S.GT.TOT) GO TO 70   
         IER = 601
         E2(1) = K
         IF (LEVEL.GE.1) WRITE (NOUT,150) K     
  150    FORMAT ('0','*** W A R N I N G ************'/'0',&
     &      '    IN ITPACK ROUTINE EQRT1S  '/' ', &
     &      '    SUCCESSIVE ITERATES TO THE',I10/' ',     &
     &      '    EIGENVALUE WERE NOT MONOTONE INCREASING ') 
!       
!                                  SET ERROR -- IRREGULAR END       
!                                  DEFLATE MINIMUM DIAGONAL ELEMENT 
!       
         S = 0.0D0
         DELTA = QP 
         DO 160 J = K,N     
            IF (D(J).GT.DELTA) GO TO 160
            I = J 
            DELTA = D(J)    
  160    CONTINUE 
!       
!                                  CONVERGENCE  
!       
  170    IF (I.LT.N) E2(I+1) = E2(I)*F/QP       
         IF (I.EQ.K) GO TO 190
         K1 = I-K 
         DO 180 JJ = 1,K1   
            J = I-JJ
            D(J+1) = D(J)-S 
            E2(J+1) = E2(J) 
  180    CONTINUE 
  190    D(K) = TOT 
         ERR = ERR+DABS(DELTA)
         E2(K) = ERR
  200 CONTINUE    
      IF (IER.EQ.0) GO TO 220 
  210 CONTINUE    
  220 IERR = IER  
      RETURN      
      END 
      INTEGER FUNCTION IPSTR (OMEGA)  
!       
!     FINDS THE SMALLEST INTEGER, IPSTR, GREATER THAN 5 SUCH THAT   
!          IPSTR * (OMEGA-1)**(IPSTR-1) .LE. 0.50. IPSTR WILL BE SET
!          IN LOOP. 
!       
! ... PARAMETER LIST:       
!       
!          OMEGA  RELAXATION FACTOR FOR SOR METHOD
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      DOUBLE PRECISION OMEGA
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER IP  
      DOUBLE PRECISION WM1  
!       
      WM1 = OMEGA-1.D0      
!       
      DO 10 IP = 6,940      
         IF (DBLE(FLOAT(IP))*(WM1**(IP-1)).GT.0.50D0) GO TO 10      
         IPSTR = IP 
         RETURN   
   10 CONTINUE    
      IPSTR = 940 
      RETURN      
!       
      END











      SUBROUTINE ITERM (NN,A,U,WK,IMTHDD)       
!       
!     THIS ROUTINE PRODUCES THE ITERATION SUMMARY LINE AT THE END   
!     OF EACH ITERATION. IF LEVEL = 5, THE LATEST APPROXIMATION     
!     TO THE SOLUTION WILL BE PRINTED.
!       
! ... PARAMETER LIST:       
!       
!          NN     ORDER OF SYSTEM OR, FOR REDUCED SYSTEM  
!                    ROUTINES, ORDER OF BLACK SUBSYSTEM   
!          A      ITERATION MATRIX    
!          U      SOLUTION ESTIMATE   
!          WK     WORK ARRAY OF LENGTH NN       
!          IMTHD  INDICATOR OF METHOD (=IMTHDD) 
!                    IMTHD = 1,  JCG  
!                    IMTHD = 2,  JSI  
!                    IMTHD = 3,  SOR  
!                    IMTHD = 4,  SSORCG 
!                    IMTHD = 5,  SSORSI 
!                    IMTHD = 6,  RSCG 
!                    IMTHD = 7,  RSSI 
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER NN,IMTHD      
      DOUBLE PRECISION A(1),U(NN),WK(NN)
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER I,IMTHDD,IP,N 
      DOUBLE PRECISION QTFF 
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
      N = NN      
      IMTHD = IMTHDD
!       
! ... PRINT VARIOUS PARAMETERS AFTER EACH ITERATION       
!       
      IF (LEVEL.LT.2) RETURN
      GO TO (10,110,170,210,50,10,110), IMTHD   
   10 IF (IN.GT.0) GO TO 30 
!       
! ... PRINT HEADER FOR JCG AND RSCG   
!       
      WRITE (NOUT,20)       
   20 FORMAT (////15X,'INTERMEDIATE OUTPUT AFTER EACH ITERATION'//  &
     &   ' NUMBER OF',5X,'CONVERGENCE',7X,'CME ',11X,'RHO',12X,'GAMMA'/ &
     &   ' ITERATIONS',4X,'TEST '//)  
!       
! ... PRINT SUMMARY LINE    
!       
   30 WRITE (NOUT,40) IN,STPTST,CME,RHO,GAMMA   
   40 FORMAT (4X,I5,3X,4D15.7)
      IF (LEVEL.GE.4) GO TO 250       
!       
      RETURN      
!       
   50 IF (IN.GT.0) GO TO 70 
!       
! ... PRINT HEADER FOR SSOR-SI
!       
      WRITE (NOUT,60)       
   60 FORMAT (////15X,'INTERMEDIATE OUTPUT AFTER EACH ITERATION'//  &
     &   ' NUMBER OF',4X,'CONVERGENCE',7X,'PARAMETER CHANGE TEST',10X,&
     &   'RHO',12X,'GAMMA'/' ITERATIONS',3X,'TEST ',11X,'LHS(QA)',7X, &
     &   'RHS(QT**FF)'//)   
!       
! ... PRINT SUMMARY LINE    
!       
   70 IP = IN-IS  
      IF (IMTHD.EQ.7) IP = 2*IP       
      IF (IP.LT.3) GO TO 90 
      QTFF = QT**FF 
      WRITE (NOUT,80) IN,STPTST,QA,QTFF,RHO,GAMMA 
   80 FORMAT (4X,I5,3X,5D15.7)
      IF (LEVEL.GE.4) GO TO 250       
      RETURN      
!       
   90 WRITE (NOUT,100) IN,STPTST,RHO,GAMMA      
  100 FORMAT (4X,I5,3X,D15.7,30X,2D15.7)
      IF (LEVEL.GE.4) GO TO 250       
      RETURN      
!       
  110 IF (IN.GT.0) GO TO 130
!       
! ... PRINT HEADER FOR J-SI AND RS-SI 
!       
      WRITE (NOUT,120)      
  120 FORMAT (////15X,'INTERMEDIATE OUTPUT AFTER EACH ITERATION'//  &
     &   ' NUMBER OF',4X,'CONVERGENCE',7X,'PARAMETER CHANGE TEST',10X,&
     &   'RHO'/' ITERATIONS',3X,'TEST ',11X,'LHS(QA)',7X,'RHS(QT**FF)'//&
     &   )
!       
! ... PRINT SUMMARY LINE    
!       
  130 IP = IN-IS  
      IF (IMTHD.EQ.7) IP = 2*IP       
      IF (IP.LT.3) GO TO 150
      QTFF = QT**FF 
      WRITE (NOUT,140) IN,STPTST,QA,QTFF,RHO    
  140 FORMAT (4X,I5,3X,5D15.7)
      IF (LEVEL.GE.4) GO TO 250       
      RETURN      
!       
  150 WRITE (NOUT,160) IN,STPTST,RHO  
  160 FORMAT (4X,I5,3X,D15.7,30X,D15.7) 
      IF (LEVEL.GE.4) GO TO 250       
      RETURN      
!       
! ... PRINT VARIOUS PARAMETERS AFTER EACH ITERATION FOR SOR.
!       
  170 IF (IN.GT.0) GO TO 190
!       
! ... PRINT HEADER FOR SOR  
!       
      WRITE (NOUT,180)      
  180 FORMAT (////15X,'INTERMEDIATE OUTPUT AFTER EACH ITERATION'//  &
     &   ' NUMBER OF',4X,'CONVERGENCE',6X,'CME ',9X,'OMEGA',7X,     &
     &   'SPECTRAL'/' ITERATIONS',3X,'TEST',38X,'RADIUS'//) 
!       
! ... PRINT SUMMARY LINE FOR SOR      
!       
  190 CONTINUE    
      WRITE (NOUT,200) IN,STPTST,CME,OMEGA,SPECR
  200 FORMAT (4X,I5,3X,4D14.7)
      IF (LEVEL.GE.4) GO TO 250       
!       
      RETURN      
!       
! ... PRINT VARIOUS PARAMETERS AFTER EACH ITERATION FOR SSOR-CG.    
!       
  210 IF (IN.GT.0) GO TO 230
!       
! ... PRINT HEADER FOR SSOR-CG
!       
      WRITE (NOUT,220)      
  220 FORMAT (////15X,'INTERMEDIATE OUTPUT AFTER EACH ITERATION'//  &
     &   ' NUMBER OF',4X,'CONVERGENCE',3X,' SPECTRAL',6X,'S-PRIME',9X,&
     &   'RHO',10X,'GAMMA'/' ITERATIONS',3X,'TEST ',10X,'RADIUS'//) 
!       
! ... PRINT SUMMARY LINE FOR SSOR-CG  
!       
  230 CONTINUE    
      WRITE (NOUT,240) IN,STPTST,SPECR,SPR,RHO,GAMMA      
  240 FORMAT (4X,I5,3X,5D14.7)
      IF (LEVEL.GE.4) GO TO 250       
      RETURN      
!       
  250 IF (IMTHD.GT.5) GO TO 270       
      WRITE (NOUT,260) IN   
  260 FORMAT ('0',2X,'ESTIMATE OF SOLUTION AT ITERATION ',I5)       
      GO TO 290   
  270 WRITE (NOUT,280) IN   
  280 FORMAT ('0',2X,'ESTIMATE OF SOLUTION AT BLACK POINTS ',       &
     &   'AT ITERATION ',I5)
  290 DO 300 I = 1,N
         WK(I) = U(I)/A(I)  
  300 CONTINUE    
      WRITE (NOUT,310) (WK(I),I=1,N)  
  310 FORMAT (2X,5(2X,D20.13))
      WRITE (NOUT,320)      
  320 FORMAT (//) 
!       
      RETURN      
      END











      SUBROUTINE IVFILL (N,IV,IVAL)   
!       
!     FILLS AN INTEGER VECTOR, IV, WITH AN INTEGER VALUE, IVAL.     
!       
! ... PARAMETER LIST:       
!       
!          N      INTEGER LENGTH OF VECTOR IV   
!          IV     INTEGER VECTOR      
!          IVAL   INTEGER CONSTANT THAT FILLS FIRST N LOCATIONS OF IV 
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER N,IVAL,IV(N)  
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER I,M,MP1       
!       
      IF (N.LE.0) RETURN    
!       
!     CLEAN UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 10  
!       
      M = MOD(N,10) 
      IF (M.EQ.0) GO TO 20  
      DO 10 I = 1,M 
         IV(I) = IVAL       
   10 CONTINUE    
      IF (N.LT.10) RETURN   
!       
   20 MP1 = M+1   
      DO 30 I = MP1,N,10    
         IV(I) = IVAL       
         IV(I+1) = IVAL     
         IV(I+2) = IVAL     
         IV(I+3) = IVAL     
         IV(I+4) = IVAL     
         IV(I+5) = IVAL     
         IV(I+6) = IVAL     
         IV(I+7) = IVAL     
         IV(I+8) = IVAL     
         IV(I+9) = IVAL     
   30 CONTINUE    
!       
      RETURN      
      END 
      SUBROUTINE OMEG (DNRM,IFLAG)    
!       
!     COMPUTES NEW VALUES FOR  CME, OMEGA, AND SPECR FOR  
!     FULLY ADAPTIVE SSOR METHODS.    
!       
! ... PARAMETER LIST:       
!       
!          DNRM   NUMERATOR OF RAYLEIGH QUOTIENT
!          IFLAG  INDICATOR OF APPROPRIATE ENTRY POINT    
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IFLAG 
      DOUBLE PRECISION DNRM 
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      DOUBLE PRECISION TEMP,ZM1,ZM2   
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
!       
      ZM1 = 0.D0  
      ZM2 = 0.D0  
      IF (IFLAG.EQ.1) GO TO 10
!       
! ... IFLAG .NE. 1, COMPUTE NEW ESTIMATE FOR CME
!       
      ZM1 = ((1.D0-SPR)*(1.D0+BETAB*OMEGA**2)-OMEGA*(2.D0-OMEGA))/(OMEGA&
     &   *(OMEGA-1.D0-SPR)) 
!       
      IF (.NOT.CASEII) ZM2 = DNRM/BDELNM
      IF (CASEII) ZM2 = DSQRT(DABS(DNRM/BDELNM))
      CME = DMAX1(CME,ZM1,ZM2)
!       
! ... IFLAG = 1, OR CONTINUATION OF IFLAG .NE. 1
!       
!        COMPUTE NEW VALUES OF OMEGA AND SPECR BASED ON CME AND BETAB 
!       
   10 IS = IN+1   
      DELSNM = DELNNM       
      IF (CME.GE.(4.D0*BETAB)) GO TO 30 
!       
! ... CME .LT. 4.D0*BETAB   
!       
      TEMP = DSQRT(DABS(1.D0-2.D0*CME+4.D0*BETAB))
      OMEGA = DMAX1((2.D0/(1.D0+TEMP)),1.D0)    
      TEMP = (1.D0-CME)/TEMP
      SPECR = (1.D0-TEMP)/(1.D0+TEMP) 
      IF (DABS(OMEGA-1.D0).LT.DRELPR) SPECR = 0.D0
      IF (LEVEL.GE.2) WRITE (NOUT,20) IN,BETAB,ZM1,ZM2,CME,OMEGA,SPECR
   20 FORMAT (/30X,'PARAMETERS WERE CHANGED AT ITERATION NO.',I5/35X, &
     &   'NEW ESTIMATE OF BETAB            =',D15.7/35X,  &
     &   'SOLUTION TO CHEBYSHEV EQN.       =',D15.7/35X,  &
     &   'SOLUTION TO RAYLEIGH QUOTIENT    =',D15.7/35X,  &
     &   'NEW ESTIMATE FOR CME             =',D15.7/35X,  &
     &   'NEW ESTIMATE FOR OMEGA           =',D15.7/35X,  &
     &   'NEW ESTIMATE FOR SPECTRAL RADIUS =',D15.7/)     
!       
      RETURN      
!       
! ... CME .GE. 4.D0*BETAB   
!       
! ... OMEGA-STAR WILL BE CHOSEN       
!       
   30 CME = 2.D0*DSQRT(DABS(BETAB))   
      OMEGA = 2.D0/(1.D0+DSQRT(DABS(1.D0-4.D0*BETAB)))    
      SPECR = OMEGA-1.D0    
      ADAPT = .FALSE.       
      PARTAD = .FALSE.      
      IF (LEVEL.GE.2) WRITE (NOUT,20) IN,BETAB,ZM1,ZM2,CME,OMEGA,SPECR
!       
      RETURN      
      END 
      LOGICAL FUNCTION OMGCHG (NDUMMY)
!       
! ... THIS FUNCTION TESTS TO SEE WHETHER OMEGA SHOULD BE CHANGED    
! ... FOR SSOR CG METHOD.   
!       
! ... PARAMETER LIST:       
!       
!          NDUMMY ARBITRARY INTEGER PARAMETER   
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER NDUMMY
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      DOUBLE PRECISION DEL1,DEL2,X    
!       
      DOUBLE PRECISION PHI  
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
!       
! ... STATEMENT FUNCTION PHI(X)       
!       
      PHI(X) = (1.D0-DSQRT(DABS(1.D0-X)))/(1.D0+DSQRT(DABS(1.D0-X)))
!       
      OMGCHG = .FALSE.      
      IF (IN-IS.LT.3) RETURN
      IF (SPECR.EQ.0.D0) GO TO 10     
      IF (SPECR.GE.SPR) RETURN
      DEL1 = -DLOG(DABS(PHI(SPECR)/PHI(SPECR/SPR)))       
      DEL2 = -DLOG(DABS(PHI(SPR)))    
      IF ((DEL1/DEL2).GE.FF) RETURN   
!       
   10 OMGCHG = .TRUE.       
!       
      RETURN      
      END







      LOGICAL FUNCTION OMGSTR (NDUMMY)
!       
!     TESTS FOR FULLY ADAPTIVE SSOR METHODS WHETHER OMEGA-STAR      
!     SHOULD BE USED FOR OMEGA AND THE ADAPTIVE PROCESS TURNED      
!     OFF.
!       
! ... PARAMETER LIST:       
!       
!          NDUMMY ARBITRARY INTEGER PARAMETER   
!       
! ... SPECIFICATION FOR ARGUMENT      
!       
      INTEGER NDUMMY
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      DOUBLE PRECISION OMSTAR,TEMP,TEMP1,X      
!       
      DOUBLE PRECISION PHI  
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
!       
! ... STATEMENT FUNCTION PHI(X)       
!       
      PHI(X) = (1.D0-DSQRT(DABS(1.D0-X)))/(1.D0+DSQRT(DABS(1.D0-X)))
!       
      OMGSTR = .FALSE.      
      IF (BETAB.GE..25D0.OR..NOT.ADAPT) RETURN  
      OMSTAR = 2.D0/(1.D0+DSQRT(DABS(1.D0-4.D0*BETAB)))   
!       
! ... TEST TO CHOSE OMEGA-STAR
!       
      IF ((OMSTAR.LE.1.D0).OR.(SPECR.LE.0.D0)) GO TO 10   
      TEMP = DLOG(DABS(PHI(OMSTAR-1.D0)))       
      TEMP1 = DLOG(DABS(PHI(SPECR)))  
      IF ((TEMP/TEMP1).LT.FF) RETURN  
!       
! ... OMEGA-STAR WAS CHOSEN 
!       
   10 OMEGA = OMSTAR
      SPECR = OMEGA-1.D0    
      OMGSTR = .TRUE.       
      ADAPT = .FALSE.       
      PARTAD = .FALSE.      
      CME = 2.D0*DSQRT(DABS(BETAB))   
      RRR = PHI(1.D0-SPECR)**2
      GAMMA = 2.D0/(2.D0-SPECR)       
      SIGE = SPECR/(2.D0-SPECR)       
      RHO = 1.D0  
      IS = IN+1   
      DELSNM = DELNNM       
      IF (LEVEL.GE.2) WRITE (NOUT,20) IN,CME,OMEGA,SPECR  
   20 FORMAT (/30X,'OMEGA-STAR, AN ALTERNATE ESTIMATE OF',&
     &   ' OMEGA, WAS CHOSEN AT ITERATION',I5/35X,&
     &   'NEW ESTIMATE FOR CME             =',D15.7/35X,  &
     &   'NEW ESTIMATE FOR OMEGA           =',D15.7/35X,  &
     &   'NEW ESTIMATE FOR SPECTRAL RADIUS =',D15.7/)     
!       
      RETURN      
      END









      SUBROUTINE PARCON (DTNRM,C1,C2,C3,C4,GAMOLD,RHOTMP,IBMTH)     
!       
!     COMPUTES ACCELERATION PARAMETERS FOR CONJUGATE GRADIENT       
!     ACCELERATED METHODS.  
!       
! ... PARAMETER LIST:       
!       
!          DTNRM  INNER PRODUCT OF RESIDUALS    
!          C1     OUTPUT: RHO*GAMMA   
!          C2     OUTPUT: RHO 
!          C3     OUTPUT: 1-RHO       
!          C4     OUTPUT: RHO*(1-GAMMA) 
!          GAMOLD OUTPUT: VALUE OF GAMMA AT PRECEDING ITERATION     
!          RHOTMP LAST ESTIMATE FOR VALUE OF RHO
!          IBMTH  INDICATOR OF BASIC METHOD BEING ACCELERATED BY CG 
!                      IBMTH = 1,   JACOBI      
!                            = 2,   REDUCED SYSTEM
!                            = 3,   SSOR
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IBMTH 
      DOUBLE PRECISION DTNRM,C1,C2,C3,C4,GAMOLD,RHOTMP    
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER IP  
      DOUBLE PRECISION RHOOLD 
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
!       
      IP = IN-IS  
!       
! ... SET RHOOLD AND GAMOLD 
!       
      RHOOLD = RHO
      GAMOLD = GAMMA
!       
! ... COMPUTE GAMMA (IN+1)  
!       
! ... FOR JACOBI OR REDUCED SYSTEM CG 
!       
      IF (IBMTH.LE.2) GAMMA = 1.D0/(1.D0-DTNRM/DELNNM)    
!       
! ... FOR SSOR CG 
!       
      IF (IBMTH.EQ.3) GAMMA = DELNNM/DTNRM      
!       
! ... COMPUTE RHO (IN+1)    
!       
      RHO = 1.D0  
      IF (IP.EQ.0) GO TO 20 
      IF (ISYM.EQ.0) GO TO 10 
      RHO = 1.D0/(1.D0-GAMMA*RHOTMP/DELSNM)     
      GO TO 20    
   10 RHO = 1.D0/(1.D0-GAMMA*DELNNM/(GAMOLD*DELSNM*RHOOLD)) 
!       
! ... COMPUTE CONSTANTS C1, C2, C3, AND C4      
!       
   20 DELSNM = DELNNM       
      RHOTMP = RHOOLD       
      C1 = RHO*GAMMA
      C2 = RHO    
      C3 = 1.D0-RHO 
      C4 = RHO*(1.D0-GAMMA) 
!       
      RETURN      
      END











      SUBROUTINE PARSI (C1,C2,C3,IBMTH) 
!       
!     COMPUTES ACCELERATION PARAMETERS FOR SEMI-ITERATIVE 
!     ACCELERATED METHODS.  
!       
! ... PARAMETER LIST:       
!       
!          C1,C2  
!           AND   
!           C3    OUTPUT ACCELERATION PARAMETERS
!          IBMTH  INDICATOR OF BASIC METHOD BEING ACCELERATED BY SI 
!                      IBMTH = 1, JACOBI
!                            = 2, REDUCED SYSTEM
!                            = 3, SSOR
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IBMTH 
      DOUBLE PRECISION C1,C2,C3       
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER IP  
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
!       
      IP = IN-IS  
      IF (IP.EQ.0) GO TO 30 
      IF (IP.EQ.1) GO TO 10 
      RHO = 1.D0/(1.D0-SIGE*SIGE*RHO*.25D0)     
      GO TO 20    
   10 RHO = 1.D0/(1.D0-SIGE*SIGE*.5D0)
!       
   20 C1 = RHO*GAMMA
      C2 = RHO    
      C3 = 1.D0-RHO 
!       
      RETURN      
!       
! ... NONADAPTIVE INITIALIZATION FOR SEMI-ITERATIVE METHODS 
!       
   30 CONTINUE    
      GO TO (40,50,60), IBMTH 
!       
! ... JSI 
!       
   40 IF (CASEII) SME = -CME
      GAMMA = 2.D0/(2.D0-CME-SME)     
      SIGE = (CME-SME)/(2.D0-CME-SME) 
      GO TO 70    
!       
! ... REDUCED SYSTEM SI     
!       
   50 GAMMA = 2.D0/(2.D0-CME*CME)     
      SIGE = CME*CME/(2.D0-CME*CME)   
      RRR = (1.D0-DSQRT(DABS(1.D0-CME*CME)))/(1.D0+DSQRT(DABS(1.D0-CME* &
     &   CME)))   
      GO TO 70    
!       
! ... SSORSI      
!       
   60 GAMMA = 2.D0/(2.D0-SPECR)       
      SIGE = SPECR/(2.D0-SPECR)       
      RRR = (1.D0-DSQRT(DABS(1.D0-SIGE*SIGE)))/(1.D0+DSQRT(DABS(1.D0- &
     &   SIGE*SIGE)))       
!       
   70 RHO = 1.D0  
      C1 = GAMMA  
      C2 = 1.D0   
      C3 = 0.D0   
!       
      RETURN      
      END

      DOUBLE PRECISION FUNCTION PBETA (NN,IA,JA,A,V,W1,W2)
!       
!     ... COMPUTES THE NUMERATOR FOR THE COMPUTATION OF BETAB IN    
!     ...  SSOR METHODS.    
!       
! ... PARAMETER LIST:       
!       
!          N      DIMENSION OF MATRIX (= NN)    
!          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
!          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
!          W1,W2  WORKSPACE VECTORS OF LENGTH N 
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(1),JA(1),NN
      DOUBLE PRECISION A(1),V(NN),W1(NN),W2(NN) 
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER I,IBGN,IEND,II,ITMP,JAI,JAJJ,JJ,K,N,NM1     
      DOUBLE PRECISION SUM,TEMP1,TEMP2
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
!       
      N = NN      
      PBETA = 0.D0
      IF (ISYM.EQ.0) GO TO 110
!       
!     ************** NON - SYMMETRIC SECTION ********************   
!       
      DO 10 I = 1,N 
         W1(I) = V(I)       
   10 CONTINUE    
      TEMP1 = 0.D0
      TEMP2 = 0.D0
      ITMP = 2    
      IBGN = IA(1)
      IEND = IA(ITMP)-1     
      IF (IEND.LT.IBGN) GO TO 30      
      DO 20 I = IBGN,IEND   
         JAI = JA(I)
         TEMP1 = TEMP1-A(I)*W1(JAI)   
   20 CONTINUE    
   30 W1(1) = TEMP1 
      W2(1) = 0.D0
      NM1 = N-1   
      DO 70 K = 2,NM1       
         TEMP1 = 0.D0       
         TEMP2 = 0.D0       
         IBGN = IA(K)       
         IEND = IA(K+1)-1   
         IF (IEND.LT.IBGN) GO TO 60   
         DO 50 I = IBGN,IEND
            JAI = JA(I)     
            IF (JAI.GT.K) GO TO 40    
            TEMP2 = TEMP2-A(I)*W1(JAI)
            GO TO 50
   40       TEMP1 = TEMP1-A(I)*W1(JAI)
   50    CONTINUE 
   60    W1(K) = TEMP1      
         W2(K) = TEMP2      
   70 CONTINUE    
      TEMP2 = 0.D0
      IBGN = IA(N)
      IEND = IA(N+1)-1      
      IF (IEND.LT.IBGN) GO TO 90      
      DO 80 I = IBGN,IEND   
         JAI = JA(I)
         TEMP2 = TEMP2-A(I)*W1(JAI)   
   80 CONTINUE    
   90 W2(N) = TEMP2 
      DO 100 I = 1,N
         PBETA = PBETA+V(I)*W2(I)     
  100 CONTINUE    
      RETURN      
!       
!     **************** SYMMETRIC SECTION *************************  
!       
  110 DO 130 II = 1,N       
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = 0.D0 
         IF (IBGN.GT.IEND) GO TO 130  
         DO 120 JJ = IBGN,IEND
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*V(JAJJ)   
  120    CONTINUE 
         PBETA = PBETA+SUM*SUM
  130 CONTINUE    
      RETURN      
!       
      END








      SUBROUTINE PBSOR (NN,IA,JA,A,U,RHS)       
!       
!     ... THIS SUBROUTINE COMPUTES A BACKWARD SOR SWEEP.  
!       
! ... PARAMETER LIST:       
!       
!          N      ORDER OF SYSTEM (= NN)
!          OMEGA  RELAXATION FACTOR   
!          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
!          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
!          U      LATEST ESTIMATE OF SOLUTION   
!          RHS    RIGHT HAND SIDE OF MATRIX PROBLEM       
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(1),JA(1),NN
      DOUBLE PRECISION A(1),U(NN),RHS(NN)       
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER I,IBGN,IEND,II,JAJJ,JJ,N,NPL1     
      DOUBLE PRECISION OMM1,SUM,UI    
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
!       
      N = NN      
      NPL1 = N+1  
      OMM1 = OMEGA-1.D0     
      IF (ISYM.EQ.0) GO TO 40 
!       
!     *************** NON - SYMMETRIC SECTION **********************
!       
      DO 30 I = 1,N 
         II = NPL1-I
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = RHS(II)      
         IF (IBGN.GT.IEND) GO TO 20   
         DO 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*U(JAJJ)   
   10    CONTINUE 
   20    U(II) = OMEGA*SUM-OMM1*U(II) 
   30 CONTINUE    
      RETURN      
!       
!     ***************** SYMMETRIC SECTION **************************
!       
   40 DO 60 II = 1,N
         UI = U(II) 
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IF (IBGN.GT.IEND) GO TO 60   
         DO 50 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            RHS(JAJJ) = RHS(JAJJ)-A(JJ)*UI      
   50    CONTINUE 
   60 CONTINUE    
!       
      DO 90 I = 1,N 
         II = NPL1-I
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = RHS(II)      
         IF (IBGN.GT.IEND) GO TO 80   
         DO 70 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*U(JAJJ)   
   70    CONTINUE 
   80    U(II) = OMEGA*SUM-OMM1*U(II) 
   90 CONTINUE    
      RETURN      
!       
      END











      SUBROUTINE PERMAT (NN,IA,JA,A,P,NEWIA,ISYM,LEVEL,NOUT,IERR)   
!       
!*********************************************************************
!       
! ... SUBROUTINE PERMAT TAKES THE SPARSE MATRIX REPRESENTATION      
!     OF THE MATRIX STORED IN THE ARRAYS IA, JA, AND A AND
!     PERMUTES BOTH ROWS AND COLUMNS OVERWRITING THE PREVIOUS       
!     STRUCTURE.  
!       
! ... PARAMETER LIST:       
!       
!         N      ORDER OF SYSTEM (= NN) 
!         IA,JA  INTEGER ARRAYS OF THE SPARSE MATRIX REPRESENTATION 
!         A      D.P. ARRAY OF THE SPARSE MATRIX REPRESENTATION     
!         P      PERMUTATION VECTOR   
!         NEWIA  INTEGER WORK VECTOR OF LENGTH N
!         ISYM   SYMMETRIC/NONSYMMETRIC STORAGE SWITCH    
!         LEVEL  SWITCH CONTROLLING LEVEL OF OUTPUT       
!         NOUT OUTPUT UNIT NUMBER     
!         IER    OUTPUT ERROR FLAG (= IERR)     
!       
!                   IER =   0  NORMAL RETURN    
!                   IER = 301  NO ENTRY IN ITH ROW OF ORIGINAL      
!                              MATRIX. IF LEVEL IS GREATER THAN     
!                              0, I WILL BE PRINTED       
!                   IER = 302  THERE IS NO ENTRY IN THE ITH ROW     
!                              OF THE PERMUTED MATRIX     
!                   IER = 303  ERROR RETURN FROM QSORT IN 
!                              SORTING THE ITH ROW OF THE 
!                              PERMUTED MATRIX  
! ... IT IS ASSUMED THAT THE I-TH ENTRY OF THE PERMUTATION VECTOR   
!     P INDICATES THE ROW THE I-TH ROW GETS MAPPED INTO.  (I.E.     
!     IF ( P(I) = J ) ROW I GETS MAPPED INTO ROW J.)      
!       
! ... THE ARRAY NEWIA IS AN INTEGER WORK VECTOR OF LENGTH N WHICH   
!     KEEPS TRACK OF WHERE THE ROWS BEGIN IN THE PERMUTED STRUCTURE.
!       
! ... PERMAT IS CAPABLE OF PERMUTING BOTH THE SYMMETRIC AND NON-    
!     SYMMETRIC FORM OF IA, JA, AND A.  IF ( ISYM .EQ. 0 ) SYMMETRIC
!     FORM IS ASSUMED.      
!       
! ... TWO EXTERNAL MODULES ARE USED BY PERMAT.  THE FIRST IS INTEGER
!     FUNCTION BISRCH WHICH USES A BISECTION SEARCH ( ORDER LOG-BASE-2
!     OF N+1 ) THROUGH THE ARRAY IA TO FIND THE ROW INDEX OF AN ARBI- 
!     TRARY ENTRY EXTRACTED FROM THE ARRAY JA. THE SECOND IS SUBROUTINE 
!     QSORT WHICH PERFORMS A QUICK SORT TO PLACE THE ENTRIES IN     
!     THE PERMUTED ROWS IN COLUMN ORDER.
!       
!*********************************************************************
!       
      INTEGER NN,IA(1),JA(1),P(NN),NEWIA(NN),ISYM,IERR    
      DOUBLE PRECISION A(1) 
!       
! ... INTERNAL VARIABLES    
!       
      INTEGER BISRCH,I,IBGN,IEND,IP,IPP,J,JAJ,JP,IER,K,N,NELS,NEXT,NPL1 
!       
      DOUBLE PRECISION SAVE,TEMP      
!       
!*********************************************************************
!       
! ... PREPROCESSING PHASE   
!       
! ...... DETERMINE THE NUMBER OF NONZEROES IN THE ROWS OF THE PERMUTED
!        MATRIX AND STORE THAT IN NEWIA.  THEN SWEEP THRU NEWIA TO MAKE 
!        NEWIA(I) POINT TO THE BEGINNING OF EACH ROW IN THE PERMUTED
!        DATA STRUCTURE.  ALSO NEGATE ALL THE ENTRIES IN JA TO INDICATE 
!        THAT THOSE ENTRIES HAVE NOT BEEN MOVED YET.      
!       
      N = NN      
      IER = 0     
      NPL1 = N+1  
      NELS = IA(NPL1)-1     
      DO 10 I = 1,N 
         NEWIA(I) = 0       
   10 CONTINUE    
      DO 30 I = 1,N 
         IP = P(I)
         IBGN = IA(I)       
         IEND = IA(I+1)-1   
         IF (IBGN.GT.IEND) GO TO 90   
         DO 20 J = IBGN,IEND
            IPP = IP
            JAJ = JA(J)     
            JP = P(JAJ)     
            IF (ISYM.EQ.0.AND.IP.GT.JP) IPP = JP
            NEWIA(IPP) = NEWIA(IPP)+1 
            JA(J) = -JAJ    
   20    CONTINUE 
   30 CONTINUE    
      IBGN = 1    
      DO 40 I = 1,N 
         K = IBGN+NEWIA(I)  
         NEWIA(I) = IBGN    
         IBGN = K 
   40 CONTINUE    
!       
! ...... PREPROCESSING NOW FINISHED.  
!       
! ...... NOW PERMUTE JA AND A.  THIS PERMUTATION WILL PERFORM THE   
!        FOLLOWING STEPS    
!       
!           1.  FIND THE FIRST ENTRY IN JA NOT PERMUTED WHICH IS    
!               INDICATED BY AN NEGATIVE VALUE IN JA      
!           2.  COMPUTE WHICH ROW THE CURRENT ENTRY IS IN.  THIS    
!               IS COMPUTED BY A BISECTION SEARCH THRU THE ARRAY    
!               IA. 
!           3.  USING THE PERMUTATION ARRAY P AND THE ARRAY NEWIA   
!               COMPUTE WHERE THE CURRENT ENTRY IS TO BE PLACED.    
!           4.  THEN PICK UP THE ENTRY WHERE THE CURRENT ENTRY WILL 
!               GO.  PUT THE CURRENT ENTRY IN PLACE.  THEN MAKE THE 
!               DISPLACED ENTRY THE CURRENT ENTRY AND LOOP TO STEP 2. 
!           5.  THIS PROCESS WILL END WHEN THE NEXT ENTRY HAS ALREADY 
!               BEEN MOVED.  THEN LOOP TO STEP 1. 
!       
      DO 70 J = 1,NELS      
         IF (JA(J).GT.0) GO TO 70     
         JAJ = -JA(J)       
         SAVE = A(J)
         NEXT = J 
         JA(J) = JAJ
!       
   50    JP = P(JAJ)
         I = BISRCH(NPL1,IA,NEXT)     
         IP = P(I)
         IPP = IP 
         IF (ISYM.NE.0.OR.IP.LE.JP) GO TO 60    
         IPP = JP 
         JP = IP  
   60    NEXT = NEWIA(IPP)  
!       
         TEMP = SAVE
         SAVE = A(NEXT)     
         A(NEXT) = TEMP     
!       
         JAJ = -JA(NEXT)    
         JA(NEXT) = JP      
         NEWIA(IPP) = NEWIA(IPP)+1    
         IF (JAJ.GT.0) GO TO 50       
!       
   70 CONTINUE    
!       
! ...... THE MATRIX IS NOW PERMUTED BUT THE ROWS MAY NOT BE IN      
!        ORDER.  THE REMAINDER OF THIS SUBROUTINE PERFORMS
!        A QUICK SORT ON EACH ROW TO SORT THE ENTRIES IN  
!        COLUMN ORDER.  THE IA ARRAY IS ALSO CORRECTED FROM 
!        INFORMATION STORED IN THE NEWIA ARRAY.  NEWIA(I) NOW       
!        POINTS TO THE FIRST ENTRY OF ROW I+1.  
!       
      IA(1) = 1   
      DO 80 I = 1,N 
         IA(I+1) = NEWIA(I) 
         K = IA(I+1)-IA(I)  
         IF (K.EQ.1) GO TO 80 
         IF (K.LT.1) GO TO 110
!       
         IBGN = IA(I)       
         CALL QSORT (K,JA(IBGN),A(IBGN),IER)    
         IF (IER.NE.0) GO TO 130      
!       
   80 CONTINUE    
!       
! ...... END OF MATRIX PERMUTATION    
!       
      GO TO 150   
!       
! ... ERROR TRAPS 
!       
! ...... NO ENTRY IN ROW I IN THE ORIGINAL SYSTEM 
!       
   90 IER = 301   
      IF (LEVEL.GE.0) WRITE (NOUT,100) I
  100 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    IN ITPACK ROUTINE PERMAT  '/' ','    NO ENTRY IN ROW ',I10&
     &   ,' OF ORIGINAL MATRIX ')     
      GO TO 150   
!       
! ...... NO ENTRY IN ROW I IN THE PERMUTED SYSTEM 
!       
  110 IER = 302   
      IF (LEVEL.GE.0) WRITE (NOUT,120) I
  120 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    IN ITPACK ROUTINE PRBNDX  '/' ','    NO ENTRY IN ROW ',I10&
     &   ,' OF PERMUTED MATRIX ')     
      GO TO 150   
!       
! ...... ERROR RETURN FROM SUBROUTINE QSORT     
!       
  130 IER = 303   
      IF (LEVEL.GE.0) WRITE (NOUT,140) I
  140 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    IN ITPACK ROUTINE QSORT   '/' ',  &
     &   '    ERROR IN SORTING PERMUTED ROW ',I12/' ',    &
     &   '    CALLED FROM ITPACK ROUTINE PRBNDX   ')      
!       
  150 CONTINUE    
      IERR = IER  
      RETURN      
      END













      SUBROUTINE PERROR5 (NN,IA,JA,A,RHS,U,W,DIGTT1,DIGTT2,IDGTTS)   
!       
!     PERROR5 COMPUTES THE RESIDUAL, R = RHS - A*U.  THE USER
!     ALSO HAS THE OPTION OF PRINTING THE RESIDUAL AND/OR THE       
!     UNKNOWN VECTOR DEPENDING ON IDGTS.
!       
! ... PARAMETER LIST:       
!       
!          N      DIMENSION OF MATRIX (= NN)    
!          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
!          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
!          RHS    RIGHT HAND SIDE OF MATRIX PROBLEM       
!          U      LATEST ESTIMATE OF SOLUTION   
!          W      WORKSPACE VECTOR    
!          DIGIT1 OUTPUT: MEASURE OF ACCURACY OF STOPPING TEST (= DIGTT1
!          DIGIT2 OUTPUT: MEASURE OF ACCURACY OF SOLUTION (= DIGTT2)
!          IDGTS   PARAMETER CONTROLING LEVEL OF OUTPUT (= IDGTTS)  
!                    IF IDGTS < 1 OR IDGTS > 4, THEN NO OUTPUT.     
!                            = 1, THEN NUMBER OF DIGITS IS PRINTED, PRO-
!                                 VIDED LEVEL .GE. 1      
!                            = 2, THEN SOLUTION VECTOR IS PRINTED, PRO- 
!                                 VIDED LEVEL .GE. 1      
!                            = 3, THEN RESIDUAL VECTOR IS PRINTED, PRO- 
!                                 VIDED LEVEL .GE. 1      
!                            = 4, THEN BOTH VECTORS ARE PRINTED, PRO- 
!                                 VIDED LEVEL .GE. 1      
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(1),JA(1),NN,IDGTTS   
      DOUBLE PRECISION A(1),RHS(NN),U(NN),W(NN),DIGTT1,DIGTT2       
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER IDGTS,N       
      DOUBLE PRECISION BNRM,DIGIT1,DIGIT2,RNRM,TEMP       
!       
! ... SPECIFICATIONS FOR FUNCTION SUBPROGRAMS   
!       
      DOUBLE PRECISION DDOT 
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
!       
      N = NN      
      IDGTS = IDGTTS
      DIGIT1 = 0.D0 
      DIGIT2 = 0.D0 
      IF (N.LE.0) GO TO 40  
!       
      DIGIT1 = -DLOG10(DABS(DRELPR))  
      IF (STPTST.GT.0.D0) DIGIT1 = -DLOG10(DABS(STPTST))  
      BNRM = DDOT(N,RHS,1,RHS,1)      
      IF (BNRM.EQ.0.D0) GO TO 10      
      CALL PMULT (N,IA,JA,A,U,W)      
      CALL WEVMW (N,RHS,W)  
      RNRM = DDOT(N,W,1,W,1)
      TEMP = RNRM/BNRM      
      IF (TEMP.EQ.0.D0) GO TO 10      
      DIGIT2 = -DLOG10(DABS(TEMP))/2.D0 
      GO TO 20    
!       
   10 DIGIT2 = -DLOG10(DABS(DRELPR))  
!       
   20 IF ((IDGTS.LT.1).OR.(LEVEL.LE.0)) GO TO 40
      WRITE (NOUT,30) DIGIT1,DIGIT2   
   30 FORMAT (/6X,'APPROX. NO. OF DIGITS (EST. REL. ERROR) =',F5.1,2X,&
     &   '(DIGIT1)'/3X,'APPROX. NO. OF DIGITS (EST. REL. RESIDUAL) =',&
     &   F5.1,2X,'(DIGIT2)')
!       
      IF (IDGTS.LE.1.OR.IDGTS.GT.4) GO TO 40    
      IF (IDGTS.NE.3) CALL VOUT (N,U,2,NOUT)    
      IF (IDGTS.GE.3) CALL VOUT (N,W,1,NOUT)    
!       
   40 CONTINUE    
      DIGTT1 = DIGIT1       
      DIGTT2 = DIGIT2       
      RETURN      
      END











      SUBROUTINE PERVEC (N,V,P)       
!       
!     THIS SUBROUTINE PERMUTES A D.P. VECTOR AS DICTATED BY THE     
!     PERMUTATION VECTOR, P.  IF P(I) = J, THEN V(J) GETS V(I).     
!       
! ... PARAMETER LIST:       
!       
!          V      D.P. VECTOR OF LENGTH N       
!          P     INTEGER PERMUTATION VECTOR     
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER N,P(N)
      DOUBLE PRECISION V(N) 
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER II,NEXT,NOW   
      DOUBLE PRECISION SAVE,TEMP      
!       
      IF (N.LE.0) RETURN    
!       
      DO 20 II = 1,N
         IF (P(II).LT.0) GO TO 20     
!       
         NEXT = P(II)       
         SAVE = V(II)       
!       
   10    CONTINUE 
         IF (P(NEXT).LT.0) GO TO 20   
         TEMP = SAVE
         SAVE = V(NEXT)     
         V(NEXT) = TEMP     
!       
         NOW = NEXT 
         NEXT = P(NOW)      
         P(NOW) = -NEXT     
         GO TO 10 
!       
   20 CONTINUE    
!       
      DO 30 II = 1,N
         P(II) = -P(II)     
   30 CONTINUE    
!       
      RETURN      
      END














      SUBROUTINE PFSOR (NN,IA,JA,A,U,RHS)       
!       
!         THIS SUBROUTINE COMPUTES A FORWARD SOR SWEEP.   
!       
! ... PARAMETER LIST:       
!       
!         N       ORDER OF SYSTEM (= NN)
!          OMEGA  RELAXATION FACTOR   
!          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
!          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
!          U      LATEST ESTIMATE OF SOLUTION   
!          RHS    RIGHT HAND SIDE OF MATRIX PROBLEM       
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(1),JA(1),NN
      DOUBLE PRECISION A(1),U(NN),RHS(NN)       
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER IBGN,IEND,II,JAJJ,JJ,N  
      DOUBLE PRECISION OMM1,SUM,UI    
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
!       
      N = NN      
      OMM1 = OMEGA-1.D0     
      IF (ISYM.EQ.0) GO TO 40 
!       
!     *********** NON - SYMMETRIC SECTION *********************     
!       
      DO 30 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = RHS(II)      
         IF (IBGN.GT.IEND) GO TO 20   
         DO 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*U(JAJJ)   
   10    CONTINUE 
   20    UI = OMEGA*SUM-OMM1*U(II)    
         U(II) = UI 
   30 CONTINUE    
      RETURN      
!       
!     ************* SYMMETRIC SECTION *************************     
!       
   40 DO 80 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = RHS(II)      
         IF (IBGN.GT.IEND) GO TO 60   
         DO 50 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*U(JAJJ)   
   50    CONTINUE 
   60    UI = OMEGA*SUM-OMM1*U(II)    
         U(II) = UI 
         IF (IBGN.GT.IEND) GO TO 80   
         DO 70 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            RHS(JAJJ) = RHS(JAJJ)-A(JJ)*UI      
   70    CONTINUE 
   80 CONTINUE    
      RETURN      
!       
      END













      SUBROUTINE PFSOR1 (NN,IA,JA,A,U,RHS)      
!       
!         THIS SUBROUTINE COMPUTES A FORWARD SOR SWEEP ON U AND     
!         COMPUTES THE NORM OF THE PSEUDO-RESIDUAL VECTOR.
!       
! ... PARAMETER LIST:       
!       
!          N      ORDER OF SYSTEM (= NN)
!          OMEGA  RELAXATION FACTOR   
!          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
!          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
!          U      LATEST ESTIMATE OF SOLUTION   
!          RHS    RIGHT HAND SIDE OF MATRIX PROBLEM       
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(1),JA(1),NN
      DOUBLE PRECISION A(1),U(NN),RHS(NN)       
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER IBGN,IEND,II,JAJJ,JJ,N  
      DOUBLE PRECISION OMM1,SUM,SUMD,UI 
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
!       
      N = NN      
      OMM1 = OMEGA-1.D0     
      SUMD = 0.D0 
      IF (ISYM.EQ.0) GO TO 40 
!       
!     **************** NON - SYMMETRIC SECTION ******************   
!       
      DO 30 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = RHS(II)      
         IF (IBGN.GT.IEND) GO TO 20   
         DO 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*U(JAJJ)   
   10    CONTINUE 
   20    CONTINUE 
         UI = OMEGA*SUM-OMM1*U(II)    
         SUMD = SUMD+(UI-U(II))**2    
         U(II) = UI 
   30 CONTINUE    
      GO TO 90    
!       
!     *************** SYMMETRIC SECTION ************************    
!       
   40 DO 80 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = RHS(II)      
         IF (IBGN.GT.IEND) GO TO 60   
         DO 50 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*U(JAJJ)   
   50    CONTINUE 
   60    CONTINUE 
         UI = OMEGA*SUM-OMM1*U(II)    
         SUMD = SUMD+(UI-U(II))**2    
         U(II) = UI 
         IF (IBGN.GT.IEND) GO TO 80   
         DO 70 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            RHS(JAJJ) = RHS(JAJJ)-A(JJ)*UI      
   70    CONTINUE 
   80 CONTINUE    
!       
   90 DELNNM = DSQRT(SUMD)  
      RETURN      
!       
      END














      SUBROUTINE PJAC (NN,IA,JA,A,U,RHS)
!       
!     ... THIS SUBROUTINE PERFORMS ONE JACOBI ITERATION.  
!       
! ... PARAMETER LIST:       
!       
!          N      DIMENSION OF MATRIX (= NN)    
!          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
!          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
!          U      ESTIMATE OF SOLUTION OF A MATRIX PROBLEM
!          RHS    ON INPUT: CONTAINS THE RIGHT HAND SIDE OF 
!                    A MATRIX PROBLEM 
!                 ON OUTPUT: CONTAINS A*U + RHS 
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(NN+1),JA(3*NN),NN
      DOUBLE PRECISION A(3*NN),U(NN),RHS(NN)
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER IBGN,IEND,II,JAJJ,JJ,N  
      DOUBLE PRECISION RHSII,UII      
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
!       
      N = NN      
      IF (ISYM.EQ.0) GO TO 30 
!       
!     *************** NON - SYMMETRIC SECTION ****************      
!       
      DO 20 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IF (IBGN.GT.IEND) GO TO 20   
         RHSII = RHS(II)    
         DO 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            RHSII = RHSII-A(JJ)*U(JAJJ) 
   10    CONTINUE 
         RHS(II) = RHSII    
   20 CONTINUE    
      RETURN      
!       
!     ************** SYMMETRIC SECTION **********************       
!       
   30 DO 50 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IF (IBGN.GT.IEND) GO TO 50   
         RHSII = RHS(II)    
         UII = U(II)
         DO 40 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            RHSII = RHSII-A(JJ)*U(JAJJ) 
            RHS(JAJJ) = RHS(JAJJ)-A(JJ)*UII     
   40    CONTINUE 
         RHS(II) = RHSII    
   50 CONTINUE    
      RETURN      
!       
      END















      SUBROUTINE PMULT (NN,IA,JA,A,U,W) 
!       
!     ... THIS SUBROUTINE PERFORMS ONE MATRIX-VECTOR MULTIPLICATION.
!       
! ... PARAMETER LIST:       
!       
!          N      DIMENSION OF MATRIX (= NN)    
!          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
!          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
!          U      LATEST ESTIMATE OF SOLUTION   
!          W      ON RETURN W CONTAINS A*U      
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(NN+1),JA(3*NN),NN
      DOUBLE PRECISION A(3*NN),U(NN),W(NN)
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER IBGN,IEND,II,JJ,N       
      DOUBLE PRECISION SUM,UII,WII    
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
!       
      N = NN      
      IF (N.LE.0) RETURN    
      IF (ISYM.EQ.0) GO TO 40 
!       
!     *************** NON - SYMMETRIC SECTION **********************
!       
      DO 30 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = 0.0D0
         IF (IBGN.GT.IEND) GO TO 20   
         DO 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM+A(JJ)*U(JAJJ)   
   10    CONTINUE 
   20    W(II) = SUM
   30 CONTINUE    
      RETURN      
!       
!     ***************** SYMMETRIC SECTION **************************
!       
   40 CALL VFILL (N,W,0.D0) 
      DO 70 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         UII = U(II)
         WII = W(II)
         IF (IBGN.GT.IEND) GO TO 60   
         DO 50 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            WII = WII+A(JJ)*U(JAJJ)   
            W(JAJJ) = W(JAJJ)+A(JJ)*UII 
   50    CONTINUE 
   60    W(II) = WII
   70 CONTINUE    
      RETURN      
!       
      END
















      SUBROUTINE PRBNDX (NN,NBLACK,IA,JA,P,IP,LEVEL,NOUT,IER)       
!       
!**************************************************************     
!       
!     THIS SUBROUTINE COMPUTES THE RED-BLACK PERMUTATION  
!     VECTORS P ( AND ITS INVERSE IP ) IF POSSIBLE.       
!       
!     THE ALGORITHM IS TO MARK THE FIRST NODE AS RED (ARBITRARY).   
!     ALL OF ITS ADJACENT NODES ARE MARKED BLACK AND PLACED IN      
!     A STACK.  THE REMAINDER OF THE CODE PULLS THE FIRST NODE      
!     OFF THE TOP OF THE STACK AND TRIES TO TYPE ITS ADJACENT NODES.
!     THE TYPING OF THE ADJACENT POINT IS A FIVE WAY CASE STATEMENT 
!     WHICH IS WELL COMMENTED BELOW (SEE DO LOOP 100).    
!       
!     THE ARRAY P IS USED BOTH TO KEEP TRACK OF THE COLOR OF A NODE 
!     (RED NODE IS POSITIVE, BLACK IS NEGATIVE) BUT ALSO THE FATHER 
!     NODE THAT CAUSED THE COLOR MARKING OF THAT POINT.  SINCE      
!     COMPLETE INFORMATION ON THE ADJACENCY STRUCTURE IS HARD TO COME 
!     BY THIS FORMS A LINK TO ENABLE THE COLOR CHANGE OF A PARTIAL  
!     TREE WHEN A RECOVERABLE COLOR CONFLICT OCCURS.      
!       
!     THE ARRAY IP IS USED AS A STACK TO POINT TO THE SET OF NODES  
!     LEFT TO BE TYPED THAT ARE KNOWN TO BE ADJACENT TO THE CURRENT 
!     FATHER NODE.
!       
!*********************************************************************
!       
!     INPUT PARAMETERS      
!       
!        N      NUMBER OF NODES.  (INTEGER, SCALAR) (= NN)
!       
!        IA,JA  ADJACENCY STRUCTURE ARRAYS.  CAN BE EITHER THE      
!               SYMMETRIC OR NONSYMMETRIC FORM.  IT IS ASSUMED      
!               THAT FOR EVERY ROW WHERE ONLY ONE ELEMENT IS
!               STORED THAT ELEMENT CORRESPONDS TO THE DIAGONAL     
!               ENTRY.  THE DIAGONAL DOES NOT HAVE TO BE THE FIRST  
!               ENTRY STORED.  (INTEGER, ARRAYS)
!        LEVEL  SWITCH FOR PRINTING   
!        NOUT OUTPUT TAPE NUMBER      
!       
!     OUTPUT PARAMETERS     
!       
!        NBLACK NUMBER OF BLACK NODES.  NUMBER OF RED NODES IS      
!               N - NBLACK.  (INTEGER, SCALAR)  
!       
!        P, IP  PERMUTATION AND INVERSE PERMUTATION VECTORS.
!               (INTEGER, ARRAYS EACH OF LENGTH N)
!       
!        IER    ERROR FLAG. (INTEGER, SCALAR)   
!       
!               IER = 0, NORMAL RETURN.  INDEXING PERFORMED 
!                        SUCCESSFULLY 
!               IER =201, RED-BLACK INDEXING NOT POSSIBLE.
!       
!******************************************************************** 
!       
      INTEGER NN,NBLACK,IA(1),JA(1),P(NN),IP(NN),IER      
!       
      INTEGER FIRST,NEXT,LAST,I,OLD,YOUNG,IBGN,IEND,J,K,CURTYP,NXTTYP,&
     &   TYPE,NRED,N
!       
!-----------------------------------------------------------------------
!       
      N = NN      
      IER = 0     
!       
!        IF ( N .LE. 0 ) GO TO 8000   
!       
      DO 10 I = 1,N 
         P(I) = 0 
         IP(I) = 0
   10 CONTINUE    
!       
! ... HANDLE THE FIRST SET OF POINTS UNTIL SOME ADJACENT POINTS     
! ... ARE FOUND   
!       
      FIRST = 1   
!       
   20 P(FIRST) = FIRST      
      IF (IA(FIRST+1)-IA(FIRST).GT.1) GO TO 40  
!       
! ... SEARCH FOR NEXT ENTRY THAT HAS NOT BEEN MARKED      
!       
      IF (FIRST.EQ.N) GO TO 130       
      IBGN = FIRST+1
      DO 30 I = IBGN,N      
         IF (P(I).NE.0) GO TO 30      
         FIRST = I
         GO TO 20 
   30 CONTINUE    
      GO TO 130   
!       
! ... FIRST SET OF ADJACENT POINTS FOUND
!       
   40 NEXT = 1    
      LAST = 1    
      IP(1) = FIRST 
!       
! ... LOOP OVER LABELED POINTS INDICATED IN THE STACK STORED IN     
! ... THE ARRAY IP
!       
   50 K = IP(NEXT)
      CURTYP = P(K) 
      NXTTYP = -CURTYP      
      IBGN = IA(K)
      IEND = IA(K+1)-1      
      IF (IBGN.GT.IEND) GO TO 110     
      DO 100 I = IBGN,IEND  
         J = JA(I)
         TYPE = P(J)
         IF (J.EQ.K) GO TO 100
!       
!================================================================== 
!       
!     THE FOLLOWING IS A FIVE WAY CASE STATEMENT DEALING WITH THE   
!     LABELING OF THE ADJACENT NODE.  
!       
! ... CASE I.  IF THE ADJACENT NODE HAS ALREADY BEEN LABELED WITH   
!              LABEL EQUAL TO NXTTYP, THEN SKIP TO THE NEXT ADJACENT
!              NODE.
!       
         IF (TYPE.EQ.NXTTYP) GO TO 100
!       
! ... CASE II.  IF THE ADJACENT NODE HAS NOT BEEN LABELED YET LABEL 
!               IT WITH NXTTYP AND ENTER IT IN THE STACK  
!       
         IF (TYPE.NE.0) GO TO 60      
         LAST = LAST+1      
         IP(LAST) = J       
         P(J) = NXTTYP      
         GO TO 100
!       
! ... CASE III.  IF THE ADJACENT NODE HAS ALREADY BEEN LABELED WITH 
!                OPPOSITE COLOR AND THE SAME FATHER SEED, THEN THERE
!                IS AN IRRECOVERABLE COLOR CONFLICT.      
!       
   60    IF (TYPE.EQ.CURTYP) GO TO 160
!       
! ... CASE IV.  IF THE ADJACENT NODE HAS THE RIGHT COLOR AND A DIFFERENT
!               FATHER NODE, THEN CHANGE ALL NODES OF THE YOUNGEST FATHE
!               NODE TO POINT TO THE OLDEST FATHER SEED AND RETAIN THE
!               SAME COLORS.
!       
         IF (TYPE*NXTTYP.LT.1) GO TO 80 
         OLD = MIN0(IABS(TYPE),IABS(NXTTYP))    
         YOUNG = MAX0(IABS(TYPE),IABS(NXTTYP))  
         DO 70 J = YOUNG,N  
            IF (IABS(P(J)).EQ.YOUNG) P(J) = ISIGN(OLD,P(J)) 
   70    CONTINUE 
         CURTYP = P(K)      
         NXTTYP = -CURTYP   
         GO TO 100
!       
! ... CASE V.  IF THE ADJACENT NODE HAS THE WRONG COLOR AND A DIFFERENT 
!              FATHER NODE, THEN CHANGE ALL NODES OF THE YOUNGEST FATHER
!              NODE TO POINT TO THE OLDEST FATHER NODE ALONG WITH   
!              CHANGING THEIR COLORS.  SINCE UNTIL THIS TIME THE    
!              YOUNGEST FATHER NODE TREE HAS BEEN INDEPENDENT NO OTHER
!              COLOR CONFLICTS WILL ARISE FROM THIS CHANGE. 
!       
   80    OLD = MIN0(IABS(TYPE),IABS(NXTTYP))    
         YOUNG = MAX0(IABS(TYPE),IABS(NXTTYP))  
         DO 90 J = YOUNG,N  
            IF (IABS(P(J)).EQ.YOUNG) P(J) = ISIGN(OLD,-P(J))
   90    CONTINUE 
         CURTYP = P(K)      
         NXTTYP = -CURTYP   
!       
! ... END OF CASE STATEMENT 
!       
!================================================================== 
!       
  100 CONTINUE    
!       
! ... ADVANCE TO NEXT NODE IN THE STACK 
!       
  110 NEXT = NEXT+1 
      IF (NEXT.LE.LAST) GO TO 50      
!       
! ... ALL NODES IN THE STACK HAVE BEEN REMOVED  
!       
! ... CHECK FOR NODES NOT LABELED.  IF ANY ARE FOUND      
! ... START THE LABELING PROCESS AGAIN AT THE FIRST       
! ... NODE FOUND THAT IS NOT LABELED. 
!       
      IBGN = FIRST+1
      DO 120 I = IBGN,N     
         IF (P(I).NE.0) GO TO 120     
         FIRST = I
         GO TO 20 
  120 CONTINUE    
!       
!===================================================================
!       
! ... ALL NODES ARE NOW TYPED EITHER RED OR BLACK 
!       
! ... GENERATE PERMUTATION VECTORS    
!       
  130 NRED = 0    
      NBLACK = 0  
      DO 150 I = 1,N
         IF (P(I).LT.0) GO TO 140     
!       
!       RED POINT 
!       
         NRED = NRED+1      
         IP(NRED) = I       
         P(I) = NRED
         GO TO 150
!       
!     BLACK POINT 
!       
  140    NBLACK = NBLACK+1  
         J = N-NBLACK+1     
         IP(J) = I
         P(I) = J 
!       
  150 CONTINUE    
!       
! ... SUCCESSFUL RED-BLACK ORDERING COMPLETED   
!       
      GO TO 180   
!       
! ........ ERROR TRAPS      
!       
! ...... N .LE. 0 
!       
!8000    IER = 200
!        GO TO 9000 
!       
! ...... TYPE CONFLICT      
!       
  160 IER = 201   
      IF (LEVEL.GE.0) WRITE (NOUT,170)
  170 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    IN ITPACK ROUTINE PRBNDX  '/' ',  &
     &   '    RED-BLACK INDEXING NOT POSSIBLE') 
!       
! ... RETURN      
!       
  180 CONTINUE    
      RETURN      
      END














      SUBROUTINE PRSBLK (NNB,NNR,IA,JA,A,UR,VB) 
!       
! ... COMPUTE A BLACK-RS SWEEP ON A RED VECTOR INTO A BLACK VECTOR  
!       
! ... PARAMETER LIST:       
!       
!         NB      NUMBER OF BLACK POINTS (= NNB)
!         NR      NUMBER OF RED POINTS (= NNR)  
!          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
!          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
!          UR     ESTIMATE OF RED SOLUTION VECTOR 
!          VB     OUTPUT: PRESENT ESTIMATE OF BLACK SOLUTION
!                    VECTOR 
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(1),JA(1),NNB,NNR     
      DOUBLE PRECISION A(1),UR(NNR),VB(NNB)     
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER I,IBGN,IEND,INR,J,JAJ,NB,NR       
      DOUBLE PRECISION SUM,URI
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
!       
      NB = NNB    
      NR = NNR    
      IF (ISYM.EQ.0) GO TO 30 
!       
!     *************** NON - SYMMETRIC SECTION **********************
!       
      DO 20 I = 1,NB
         INR = I+NR 
         IBGN = IA(INR)     
         IEND = IA(INR+1)-1 
         SUM = VB(I)
         IF (IBGN.GT.IEND) GO TO 20   
         DO 10 J = IBGN,IEND
            JAJ = JA(J)     
            SUM = SUM-A(J)*UR(JAJ)    
   10    CONTINUE 
         VB(I) = SUM
   20 CONTINUE    
      RETURN      
!       
!     ***************** SYMMETRIC SECTION **************************
!       
   30 DO 50 I = 1,NR
         IBGN = IA(I)       
         IEND = IA(I+1)-1   
         IF (IBGN.GT.IEND) GO TO 50   
         URI = UR(I)
         DO 40 J = IBGN,IEND
            JAJ = JA(J)-NR  
            VB(JAJ) = VB(JAJ)-A(J)*URI
   40    CONTINUE 
   50 CONTINUE    
!       
      RETURN      
      END
















      SUBROUTINE PRSRED (NNB,NNR,IA,JA,A,UB,VR) 
!       
! ... COMPUTES A RED-RS SWEEP ON A BLACK VECTOR INTO A RED VECTOR.  
!       
! ... PARAMETER LIST:       
!       
!         NB      NUMBER OF BLACK POINTS (= NNR)
!         NR      NUMBER OF RED POINTS (= NNB)  
!          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
!          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
!          UB     PRESENT ESTIMATE OF BLACK SOLUTION VECTOR 
!          VR     OUTPUT: PRESENT ESTIMATE OF RED SOLUTION VECTOR   
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(1),JA(1),NNB,NNR     
      DOUBLE PRECISION A(1),UB(NNB),VR(NNR)     
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER IBGN,IEND,II,JAJJ,JJ,NB,NR
      DOUBLE PRECISION SUM  
!       
      NB = NNB    
      NR = NNR    
      DO 20 II = 1,NR       
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IF (IBGN.GT.IEND) GO TO 20   
         SUM = VR(II)       
         DO 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)-NR
            SUM = SUM-A(JJ)*UB(JAJJ)  
   10    CONTINUE 
         VR(II) = SUM       
   20 CONTINUE    
!       
      RETURN      
      END


















      SUBROUTINE PSSOR1 (NN,IA,JA,A,U,RHS,FR,BR)
!       
!     ... COMPUTES COMPLETE SSOR SWEEP ON U.  U IS OVERWRITTEN      
!     ... WITH THE NEW ITERANT, FR AND BR WILL CONTAIN    
!     ... THE FORWARD AND BACKWARD RESIDUALS ON OUTPUT.   
!       
! ... PARAMETER LIST:       
!       
!         N       ORDER OF SYSTEM (= NN)
!          OMEGA  RELAXATION FACTOR   
!          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
!          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
!          U      ESTIMATE OF SOLUTION
!          RHS    RIGHT HAND SIDE OF MATRIX PROBLEM       
!          FR,BR  OUTPUT: FORWARD AND BACKWARD RESIDUALS RESPECTIVELY 
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(1),JA(1),NN
      DOUBLE PRECISION A(1),U(NN),RHS(NN),FR(NN),BR(NN)   
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER I,IBGN,IEND,II,JAJJ,JJ,N,NPL1     
      DOUBLE PRECISION OMM1,SUM,UII   
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
!       
      N = NN      
      NPL1 = N+1  
      OMM1 = OMEGA-1.D0     
      IF (ISYM.EQ.0) GO TO 40 
!       
!     *************** NON - SYMMETRIC SECTION **********************
!       
!     ... FORWARD SWEEP     
!       
      DO 30 II = 1,N
         BR(II) = U(II)     
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = RHS(II)      
         IF (IBGN.GT.IEND) GO TO 20   
         DO 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*U(JAJJ)   
   10    CONTINUE 
   20    UII = OMEGA*SUM-OMM1*U(II)   
         FR(II) = UII-U(II) 
         U(II) = UII
   30 CONTINUE    
      GO TO 90    
!       
!     ***************** SYMMETRIC SECTION **************************
!       
!     ... FORWARD SWEEP     
!       
   40 DO 80 II = 1,N
         BR(II) = U(II)     
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = RHS(II)      
         IF (IBGN.GT.IEND) GO TO 60   
         DO 50 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*U(JAJJ)   
   50    CONTINUE 
   60    UII = OMEGA*SUM-OMM1*U(II)   
         FR(II) = UII-U(II) 
         U(II) = UII
         IF (IBGN.GT.IEND) GO TO 80   
         DO 70 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            RHS(JAJJ) = RHS(JAJJ)-A(JJ)*UII     
   70    CONTINUE 
   80 CONTINUE    
!       
!     ... BACKWARD SWEEP    
!       
   90 DO 120 I = 1,N
         II = NPL1-I
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         UII = RHS(II)      
         IF (IBGN.GT.IEND) GO TO 110  
         DO 100 JJ = IBGN,IEND
            JAJJ = JA(JJ)   
            UII = UII-A(JJ)*U(JAJJ)   
  100    CONTINUE 
  110    U(II) = OMEGA*UII-OMM1*U(II) 
         BR(II) = U(II)-BR(II)
  120 CONTINUE    
!       
      RETURN      
!       
      END














      SUBROUTINE PSTOP (N,U,DNRM,CCON,IFLAG,Q1) 
!       
!     THIS SUBROUTINE PERFORMS A TEST TO SEE IF THE ITERATIVE       
!     METHOD HAS CONVERGED TO A SOLUTION INSIDE THE ERROR 
!     TOLERANCE, ZETA.      
!       
! ... PARAMETER LIST:       
!       
!          N      ORDER OF SYSTEM     
!          U      PRESENT SOLUTION ESTIMATE     
!          DNRM   INNER PRODUCT OF PSEUDO-RESIDUALS AT PRECEDING    
!                    ITERATION
!          CON    STOPPING TEST PARAMETER (= CCON)
!          IFLAG  STOPPING TEST INTEGER FLAG    
!                    IFLAG = 0,  SOR ITERATION ZERO       
!                    IFLAG = 1,  NON-RS METHOD  
!                    IFLAG = 2,  RS METHOD      
!          Q1     STOPPING TEST LOGICAL FLAG    
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER N,IFLAG       
      DOUBLE PRECISION U(N),DNRM,CCON 
      LOGICAL Q1  
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      DOUBLE PRECISION CON,TL,TR,UOLD 
!       
! ... SPECIFICATIONS FOR ARGUMENT SUBROUTINES   
!       
      DOUBLE PRECISION DDOT 
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
!       
      CON = CCON  
      HALT = .FALSE.
!       
!     SPECIAL PROCEDURE FOR ZEROTH ITERATION    
!       
      IF (IN.GE.1) GO TO 10 
      Q1 = .FALSE.
      UDNM = 1.D0 
      STPTST = 1.D3 
      IF (IFLAG.LE.0) RETURN
!       
! ... TEST IF UDNM NEEDS TO BE RECOMPUTED       
!       
   10 CONTINUE    
      IF (Q1) GO TO 20      
      IF ((IN.GT.5).AND.(MOD(IN,5).NE.0)) GO TO 20
      UOLD = UDNM 
      UDNM = DDOT(N,U,1,U,1)
      IF (UDNM.EQ.0.D0) UDNM = 1.D0   
      IF ((IN.GT.5).AND.(DABS(UDNM-UOLD).LE.UDNM*ZETA)) Q1 = .TRUE. 
!       
! ... COMPUTE STOPPING TEST 
!       
   20 TR = DSQRT(UDNM)      
      TL = 1.D0   
      IF (CON.EQ.1.D0) GO TO 40       
      IF (IFLAG.EQ.2) GO TO 30
      TL = DSQRT(DNRM)      
      TR = TR*(1.D0-CON)    
      GO TO 40    
   30 TL = DSQRT(2.D0*DNRM) 
      TR = TR*(1.D0-CON*CON)
   40 STPTST = TL/TR
      IF (TL.GE.TR*ZETA) RETURN       
      HALT = .TRUE. 
!       
      RETURN      
      END 
      DOUBLE PRECISION FUNCTION PVTBV (N,IA,JA,A,V)       
!       
!     THIS FUNCTION COMPUTES  (V**T)*A*V.       
!       
! ... PARAMETER LIST:       
!       
!          N      DIMENSION OF MATRIX 
!          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
!          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
!          V      D.P. VECTOR OF LENGTH N       
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(1),JA(1),N 
      DOUBLE PRECISION A(1),V(N)      
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER IBGN,IEND,II,JAJJ,JJ    
      DOUBLE PRECISION SUM,SUMR       
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
!       
      PVTBV = 0.D0
      SUM = 0.D0  
      DO 20 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IF (IBGN.GT.IEND) GO TO 20   
         SUMR = 0.D0
         DO 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUMR = SUMR-A(JJ)*V(JAJJ) 
   10    CONTINUE 
         SUM = SUM+V(II)*SUMR 
   20 CONTINUE    
!       
      IF (ISYM.EQ.0) SUM = 2.D0*SUM   
      PVTBV = SUM 
!       
      RETURN      
      END















      SUBROUTINE QSORT (NN,KEY,DATA,ERROR)      
!       
!     ==================================================================
!       
!     Q U I C K S O R T     
!       
!         IN THE STYLE OF THE CACM PAPER BY BOB SEDGEWICK, OCTOBER 1978 
!       
!     INPUT:      
!         N    -- NUMBER OF ELEMENTS TO BE SORTED (= NN)  
!         KEY  -- AN ARRAY OF LENGTH  N  CONTAINING THE VALUES      
!                 WHICH ARE TO BE SORTED
!         DATA -- A SECOND ARRAY OF LENGTH  N  CONTAINING DATA      
!                 ASSOCIATED WITH THE INDIVIDUAL KEYS.    
!       
!     OUTPUT:     
!         KEY  -- WILL BE ARRANGED SO THAT VALUES ARE IN INCREASING 
!                 ORDER     
!         DATA -- REARRANGED TO CORRESPOND TO REARRANGED KEYS       
!         ERROR -- WILL BE ZERO UNLESS YOUR INPUT FILE WAS OF TRULY 
!                  ENORMOUS LENGTH, IN WHICH CASE IT WILL BE EQUAL TO 1.
!       
!     ==================================================================
!       
      INTEGER NN,ERROR,KEY(NN)
      DOUBLE PRECISION DATA(NN)       
!       
!     ------------------------
!       
      INTEGER TOP,LEFT,RIGHT,I,J,TINY,V,K,IP1,JM1,LLEN,RLEN,N       
      LOGICAL DONE
      DOUBLE PRECISION D    
      INTEGER STKLEN,STACK(30)
!       
      DATA TINY,STKLEN / 9,30 /       
!       
!     -----------------------------------       
!       
!     ... PROGRAM IS A DIRECT TRANSLATION INTO FORTRAN OF SEDGEWICK^S 
!         PROGRAM 2, WHICH IS NON-RECURSIVE, IGNORES FILES OF LENGTH
!         LESS THAN 'TINY' DURING PARTITIONING, AND USES MEDIAN OF THREE
!         PARTITIONING.     
!       
      N = NN      
      IF (N.EQ.1) RETURN    
      IF (N.LE.0) GO TO 240 
!       
      ERROR = 0   
      TOP = 1     
      LEFT = 1    
      RIGHT = N   
      DONE = (N.LE.TINY)    
!       
      IF (DONE) GO TO 150   
      CALL IVFILL (STKLEN,STACK,0)    
!       
!     ===========================================================   
!     QUICKSORT -- PARTITION THE FILE UNTIL NO SUBFILE REMAINS OF   
!     LENGTH GREATER THAN 'TINY'      
!     ===========================================================   
!       
!     ... WHILE NOT DONE DO ...       
!       
   10 IF (DONE) GO TO 150   
!       
!         ... FIND MEDIAN OF LEFT, RIGHT AND MIDDLE ELEMENTS OF CURRENT 
!             SUBFILE, WHICH IS  KEY(LEFT), ..., KEY(RIGHT) 
!       
      LFRH2 = (LEFT+RIGHT)/2
      K = KEY(LFRH2)
      D = DATA(LFRH2)       
      KEY(LFRH2) = KEY(LEFT)
      DATA(LFRH2) = DATA(LEFT)
      KEY(LEFT) = K 
      DATA(LEFT) = D
!       
      IF (KEY(LEFT+1).LE.KEY(RIGHT)) GO TO 20   
      K = KEY(LEFT+1)       
      D = DATA(LEFT+1)      
      KEY(LEFT+1) = KEY(RIGHT)
      DATA(LEFT+1) = DATA(RIGHT)      
      KEY(RIGHT) = K
      DATA(RIGHT) = D       
!       
   20 IF (KEY(LEFT).LE.KEY(RIGHT)) GO TO 30     
      K = KEY(LEFT) 
      D = DATA(LEFT)
      KEY(LEFT) = KEY(RIGHT)
      DATA(LEFT) = DATA(RIGHT)
      KEY(RIGHT) = K
      DATA(RIGHT) = D       
!       
   30 IF (KEY(LEFT+1).LE.KEY(LEFT)) GO TO 40    
      K = KEY(LEFT+1)       
      D = DATA(LEFT+1)      
      KEY(LEFT+1) = KEY(LEFT) 
      DATA(LEFT+1) = DATA(LEFT)       
      KEY(LEFT) = K 
      DATA(LEFT) = D
!       
   40 V = KEY(LEFT) 
!       
!         ... V IS NOW THE MEDIAN VALUE OF THE THREE KEYS.  NOW MOVE
!             FROM THE LEFT AND RIGHT ENDS SIMULTANEOUSLY, EXCHANGING 
!             KEYS AND DATA UNTIL ALL KEYS LESS THAN  V  ARE PACKED TO
!             THE LEFT, ALL KEYS LARGER THAN  V  ARE PACKED TO THE  
!             RIGHT.
!       
      I = LEFT+1  
      J = RIGHT   
!       
!         LOOP    
!             REPEAT I = I+1 UNTIL KEY(I) >= V; 
!             REPEAT J = J-1 UNTIL KEY(J) <= V; 
!         EXIT IF J < I;    
!             << EXCHANGE KEYS I AND J >>       
!         END     
!       
   50 CONTINUE    
   60 I = I+1     
      IF (KEY(I).LT.V) GO TO 60       
!       
   70 J = J-1     
      IF (KEY(J).GT.V) GO TO 70       
!       
      IF (J.LT.I) GO TO 80  
      K = KEY(I)  
      D = DATA(I) 
      KEY(I) = KEY(J)       
      DATA(I) = DATA(J)     
      KEY(J) = K  
      DATA(J) = D 
      GO TO 50    
!       
   80 K = KEY(LEFT) 
      D = DATA(LEFT)
      KEY(LEFT) = KEY(J)    
      DATA(LEFT) = DATA(J)  
      KEY(J) = K  
      DATA(J) = D 
!       
!         ... WE HAVE NOW PARTITIONED THE FILE INTO TWO SUBFILES,   
!             ONE IS (LEFT ... J-1)  AND THE OTHER IS (I...RIGHT).  
!             PROCESS THE SMALLER NEXT.  STACK THE LARGER ONE.      
!       
      LLEN = J-LEFT 
      RLEN = RIGHT-I+1      
      IF (MAX0(LLEN,RLEN).GT.TINY) GO TO 100    
!       
!             ... BOTH SUBFILES ARE TINY, SO UNSTACK NEXT LARGER FILE 
!       
      IF (TOP.EQ.1) GO TO 90
      TOP = TOP-2 
      LEFT = STACK(TOP)     
      RIGHT = STACK(TOP+1)  
      GO TO 10    
!       
   90 DONE = .TRUE. 
!       
      GO TO 10    
!       
!             ... ELSE ONE OR BOTH SUBFILES ARE LARGE     
!       
  100 IF (MIN0(LLEN,RLEN).GT.TINY) GO TO 120    
!       
!             ... ONE SUBFILE IS SMALL, ONE LARGE.  IGNORE THE SMALL ONE
!       
      IF (LLEN.GT.RLEN) GO TO 110     
      LEFT = I    
      GO TO 10    
!       
  110 RIGHT = J-1 
!       
      GO TO 10    
!       
!         ... ELSE BOTH ARE LARGER THAN TINY.  ONE MUST BE STACKED. 
!       
  120 IF (TOP.GE.STKLEN) GO TO 240    
      IF (LLEN.GT.RLEN) GO TO 130     
      STACK(TOP) = I
      STACK(TOP+1) = RIGHT  
      RIGHT = J-1 
      GO TO 140   
!       
  130 STACK(TOP) = LEFT     
      STACK(TOP+1) = J-1    
      LEFT = I    
!       
  140 TOP = TOP+2 
!       
      GO TO 10    
!       
!     ------------------------------------------------------------  
!     INSERTION SORT THE ENTIRE FILE, WHICH CONSISTS OF A LIST      
!     OF 'TINY' SUBFILES, LOCALLY OUT OF ORDER, GLOBALLY IN ORDER.  
!     ------------------------------------------------------------  
!       
!     ... FIRST, FIND LARGEST ELEMENT IN 'KEY'  
!       
  150 I = N-1     
      LEFT = MAX0(0,N-TINY) 
      K = KEY(N)  
      J = N       
!       
  160 IF (I.LE.LEFT) GO TO 180
      IF (KEY(I).LE.K) GO TO 170      
      K = KEY(I)  
      J = I       
!       
  170 I = I-1     
      GO TO 160   
!       
  180 IF (J.EQ.N) GO TO 190 
!       
!     ... LARGEST ELEMENT WILL BE IN  KEY(N)    
!       
      KEY(J) = KEY(N)       
      KEY(N) = K  
      D = DATA(N) 
      DATA(N) = DATA(J)     
      DATA(J) = D 
!       
!     ... INSERTION SORT ... FOR I := N-1 STEP -1 TO 1 DO ...       
!       
  190 I = N-1     
      IP1 = N     
!       
  200 IF (KEY(I).LE.KEY(IP1)) GO TO 220 
!       
!             ... OUT OF ORDER ... MOVE UP TO CORRECT PLACE 
!       
      K = KEY(I)  
      D = DATA(I) 
      J = IP1     
      JM1 = I     
!       
!             ... REPEAT ... UNTIL 'CORRECT PLACE FOR K FOUND'      
!       
  210 KEY(JM1) = KEY(J)     
      DATA(JM1) = DATA(J)   
      JM1 = J     
      J = J+1     
      IF (KEY(J).LT.K) GO TO 210      
!       
      KEY(JM1) = K
      DATA(JM1) = D 
!       
  220 IP1 = I     
      I = I-1     
      IF (I.GT.0) GO TO 200 
!       
  230 RETURN      
!       
  240 ERROR = 1   
      GO TO 230   
!       
      END















      SUBROUTINE SBAGN (N,NZ,IA,JA,A,IWORK,LEVELL,NOUTT,IERR)       
!       
! ... THE ROUTINES SBINI, SBSIJ, AND SBEND CREATE A SPARSE
!     MATRIX STRUCTURE BY MEANS OF A LINKED LIST WHICH IS 
!     DESTROYED BY SBEND. SBAGN CREATES A NEW LINKED LIST 
!     SO THAT ELEMENTS MAY BE ADDED TO THE MATRIX AFTER SBEND       
!     HAS BEEN CALLED. SBAGN SHOULD BE CALLED WITH THE APPRO-       
!     PRIATE PARAMETERS, AND THEN SBSIJ AND SBEND CAN BE CALLED     
!     TO ADD THE ELEMENTS AND COMPLETE THE SPARSE MATRIX STRUC-     
!     TURE.       
!       
! ... PARAMETER LIST:       
!       
!           N       ORDER OF THE SYSTEM 
!           NZ      MAXIMUM NUMBER OF NON-ZERO ELEMENTS   
!                   IN THE SYSTEM     
!           IA, JA  INTEGER ARRAYS OF THE SPARSE
!                   MATRIX STRUCTURE  
!           A       D.P. ARRAY OF THE SPARSE MATRIX       
!                   STRUCTURE 
!           IWORK   WORK ARRAY OF DIMENSION NZ  
!           LEVEL   OUTPUT LEVEL CONTROL (= LEVELL)       
!           NOUT  OUTPUT FILE NUMBER (= NOUTT)  
!           IER     ERROR FLAG (= IERR). POSSIBLE RETURNS ARE       
!                      IER = 0, SUCCESSFUL COMPLETION     
!                          = 703, NZ TOO SMALL - NO MORE  
!                                 ELEMENTS CAN BE ADDED   
!       
! ... SPECIFICTIONS FOR ARGUMENTS     
!       
      INTEGER NZ,IA(1),JA(1),IWORK(NZ),N,LEVELL,NOUTT,IERR
      DOUBLE PRECISION A(NZ)
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER I,IER,J,LEVEL,NOUT,NADD,NADDP1,NOW,NP1,NTO,NTN
!       
! ... INITIALIZE LOCAL VARIABLES AND MAKE ERROR CHECK     
!       
      NOW = IA(N+1)-1       
      NADD = NZ-NOW 
      IER = 0     
      LEVEL = LEVELL
      NOUT = NOUTT
      IF (NADD.LE.0) IER = 703
      IF (IER.EQ.0) GO TO 20
      IF (LEVEL.GE.0) WRITE (NOUT,10) IER       
   10 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    IN ITPACK ROUTINE SBAGN   '/' ','    IER = ',I10/' ', &
     &   '    NZ TOO SMALL - NO ROOM FOR NEW ENTRY')      
      GO TO 90    
!       
! ... SHIFT ELEMENTS OF A AND JA DOWN AND ADD ZERO FILL   
!       
   20 NTO = NOW   
      NTN = NZ    
      DO 30 I = 1,NOW       
         JA(NTN) = JA(NTO)  
         A(NTN) = A(NTO)    
         NTO = NTO-1
         NTN = NTN-1
   30 CONTINUE    
      DO 40 I = 1,NADD      
         JA(I) = 0
         A(I) = 0.D0
   40 CONTINUE    
!       
! ... UPDATE IA TO REFLECT DOWNWARD SHIFT IN A AND JA     
!       
      NP1 = N+1   
      DO 50 I = 1,NP1       
         IA(I) = IA(I)+NADD 
   50 CONTINUE    
!       
! ... CREATE LINKED LIST    
!       
      NADDP1 = NADD+1       
      DO 60 I = NADDP1,NZ   
         IWORK(I) = I+1     
   60 CONTINUE    
      DO 70 I = 1,NADD      
         IWORK(I) = 0       
   70 CONTINUE    
      DO 80 I = 1,N 
         J = IA(I+1)-1      
         IWORK(J) = -I      
   80 CONTINUE    
!       
! ... INDICATE IN LAST POSITION OF IA HOW MANY SPACES     
!     ARE LEFT IN A AND JA FOR ADDITION OF ELEMENTS       
!       
      IA(N+1) = NADD
      RETURN      
!       
! ... ERROR RETURN
!       
   90 IERR = IER  
      RETURN      
      END













      SUBROUTINE SBELM (NN,IA,JA,A,RHS,IW,RW,TOL,ISYM,LEVEL,NOUT,IER) 
!       
! ... SBELM IS DESIGNED TO REMOVE ROWS AND COLUMNS OF THE MATRIX    
! ... WHERE DABS(A(I,J))/A(I,I) .LE. TOL FOR J = 1 TO N AND A(I,I)  
! ... .GT. 0. THIS IS TO TAKE CARE OF MATRICES ARISING    
! ... FROM FINITE ELEMENT DISCRETIZATIONS OF PDE^S WITH DIRICHLET   
! ... BOUNDARY CONDITIONS.  ANY SUCH ROWS AND CORRESPONDING COLUMNS 
! ... ARE THEN SET TO THE IDENTITY AFTER CORRECTING RHS.  
!       
! ... PARAMETER LIST:       
!       
!          N      DIMENSION OF MATRIX (= NN)    
!          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
!          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
!          RHS    RIGHT HAND SIDE OF MATRIX PROBLEM       
!          IW,RW  WORK ARRAYS OF LENGTH N       
!          TOL    TOLERANCE FACTOR    
!          ISYM   FLAG FOR TYPE OF STORAGE FOR SYSTEM     
!                 (0: SYMMETRIC, 1:NONSYMMETRIC)
!          LEVEL  PRINTING SWITCH FOR ERROR CONDITION     
!          NOUT OUTPUT TAPE NUMBER    
!          IER    ERROR FLAG: NONZERO VALUE ON RETURN MEANS 
!                    101 : DIAGONAL ENTRY NOT POSITIVE    
!                    102 : THERE IS NO DIAGONAL ENTRY IN ROW
!       
!********************************************************************** 
!       
!     UPDATE.  SBELM HAS BEEN REWRITTEN TO SPEED UP THE LOCATION OF 
!              OF ROWS WHICH ARE TO BE ELIMINATED.  THIS IS DONE BY 
!              FIRST STORING THE LARGEST ELEMENT OF EACH ROW IN     
!              THE ARRAY RW.  THE DIAGONAL ENTRY IS THEN COMPARED   
!              WITH THE CORRESPONDING ELEMENT IN RW.  IF IT IS      
!              DECIDED TO ELIMINATE THE ROW THEN IT IS MARKED FOR   
!              ELIMINATION. 
!       
!              WHEN A ROW IS TO BE ELIMINATED ITS DIAGONAL ENTRY    
!              IS STORED IN  RW  AND  IW IS MARKED BY A NONZERO     
!              (WHICH IS THIS ROW NUMBER)       
!       
!              ROWS WHICH HAVE ONLY DIAGONAL ENTRIES ARE NOT
!              ALTERED.     
!       
!*********************************************************************
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER NN,IA(1),JA(1),IW(NN),ISYM,LEVEL,NOUT,IER   
      DOUBLE PRECISION A(1),RHS(NN),RW(NN),TOL  
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER IBGN,ICNT,IEND,JJ,JJDI,KK,N       
      DOUBLE PRECISION DI   
!       
      N = NN      
!       
!        IF (N .GE. 1) GO TO 10       
!           IER = 100       
!           RETURN
! 10     CONTINUE 
!       
! ... STORE THE LARGEST (DABSOLUTE VALUE) OFF DIAGONAL ENTRY FOR    
! ... ROW II IN RW(II).     
!       
      IER = 0     
      ICNT = 0    
      DO 10 II = 1,N
         RW(II) = 0.0D0     
         IW(II) = 0 
   10 CONTINUE    
      DO 20 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IF (IBGN.GT.IEND) GO TO 140  
         DO 20 JJ = IBGN,IEND 
            KK = JA(JJ)     
            IF (KK.EQ.II) GO TO 20    
            RW(II) = DMAX1(RW(II),DABS(A(JJ)))  
            IF (ISYM.NE.0) GO TO 20   
            RW(KK) = DMAX1(RW(KK),DABS(A(JJ)))  
   20 CONTINUE    
!       
! ... FOR II = 1 TO N FIND THE DIAGONAL ENTRY IN ROW II   
!       
      DO 80 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         DO 40 JJ = IBGN,IEND 
            IF (JA(JJ).NE.II) GO TO 40
            DI = A(JJ)      
            JJDI = JJ       
            IF (DI.GT.0.D0) GO TO 50  
            IER = 101       
            IF (LEVEL.GE.0) WRITE (NOUT,30) II,DI 
   30       FORMAT ('0','*** F A T A L     E R R O R ************'/'0', &
     &         '    IN ITPACK ROUTINE SBELM   '/' ',      &
     &         '    DIAGONAL ELEMENT',I10,' NOT POSITIVE  '/' ',    &
     &         '    CURRENT VALUE = ',D15.8)    
            RETURN
   40    CONTINUE 
         GO TO 140
   50    CONTINUE 
!       
! ... CHECK THE SIZE OF THE LARGEST OFF DIAGONAL ELEMENT  
! ... ( STORED IN RW(II) ) AGAINST THE DIAGONAL ELEMENT DII.
!       
         IF (RW(II).NE.0.0D0) GO TO 60
         IF (1.0D0/DI.LE.TOL) GO TO 70
         GO TO 80 
   60    IF (RW(II)/DI.GT.TOL) GO TO 80 
!       
! ... THE OFF DIAGONAL ELEMENTS ARE SMALL COMPARED TO THE DIAGONAL  
! ... THEREFORE MARK IT FOR ELIMINATION AND PERFORM INITIAL 
! ... PROCESSING  
!       
   70    ICNT = ICNT+1      
         IW(II) = II
         RW(II) = DI
         A(JJDI) = 1.0D0    
         RHS(II) = RHS(II)/DI 
!       
   80 CONTINUE    
!       
! ... ELIMINATE THE ROWS AND COLUMNS INDICATED BY THE NONZERO       
! ... ENTRIES IN IW.  THERE ARE ICNT OF THEM    
!       
      IF (ICNT.EQ.0) GO TO 130
!       
! ... THE ELIMINATION IS AS FOLLOWS:  
!       
!     FOR II = 1 TO N DO    
!        IF ( IW(II) .NE. 0 ) THEN    
!           SET DIAGONAL VALUE TO 1.0  ( ALREADY DONE )   
!           SET RHS(II) = RHS(II) / RW(II)   ( ALREADY DONE )       
!           FIND NONZERO OFFDIAGONAL ENTRIES  KK
!           IF ( IW(KK) .EQ. 0 ) FIX UP RHS(KK)  WHEN USING SYMMETRIC ST
!           SET A(II,KK) = 0.0
!        ELSE ( I.E.  IW(II) .EQ. 0  )
!           FIND NONZERO OFFDIAGONAL ENTRIES   KK 
!           IF ( IW(KK) .NE. 0 ) FIX UP RHS(II) 
!                                AND SET A(II,KK) = 0.0   
!        END IF   
!     END DO      
!       
      DO 120 II = 1,N       
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IF (IW(II).EQ.0) GO TO 100   
!       
! ... THE II-TH ROW IS TO BE ELIMINATED 
!       
         DO 90 JJ = IBGN,IEND 
            KK = JA(JJ)     
            IF (KK.EQ.II) GO TO 90    
            IF ((IW(KK).EQ.0).AND.(ISYM.EQ.0)) RHS(KK) = RHS(KK)-A(JJ)* &
     &         RHS(II)      
            A(JJ) = 0.0D0   
   90    CONTINUE 
         GO TO 120
!       
! ... THE II-TH ROW IS KEPT.  CHECK THE OFF-DIAGONAL ENTRIES
!       
  100    DO 110 JJ = IBGN,IEND
            KK = JA(JJ)     
            IF (KK.EQ.II.OR.IW(KK).EQ.0) GO TO 110
            RHS(II) = RHS(II)-A(JJ)*RHS(KK)     
            A(JJ) = 0.0D0   
  110    CONTINUE 
!       
  120 CONTINUE    
!       
  130 RETURN      
!       
! ... ERROR TRAPS -- NO DIAGONAL ENTRY IN ROW II (ROW MAY BE EMPTY).
!       
  140 CONTINUE    
      IER = 102   
      IF (LEVEL.GE.0) WRITE (NOUT,150) II       
  150 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    IN ITPACK ROUTINE SBELM   '/' ',  &
     &   '    NO DIAGONAL ENTRY IN ROW  ',I10)  
!       
      RETURN      
      END












      SUBROUTINE SBEND (N,NZ,IA,JA,A,IWORK)     
!       
!***********************************************************************
!       
!     SBEND IS THE THIRD OF A SUITE OF SUBROUTINES TO AID THE       
!     USER TO CONSTRUCT THE  IA, JA, A DATA STRUCTURE USED IN       
!     ITPACK.     
!       
!     SBEND RESTRUCTURES THE LINKED LIST DATA STRUCTURE BUILT BY    
!     SBINI AND SBSIJ INTO THE FINAL DATA STRUCTURE REQUIRE BY      
!     ITPACK.  THE RESTRUCTURING CAN TAKE PLACE IN THE MINIMUM      
!     AMOUNT OF MEMORY REQUIRED TO HOLD THE NONZERO STRUCTURE OF    
!     THE SPARSE MATRIX BUT WILL RUN QUICKER IF MORE STORAGE
!     IS ALLOWED. 
!       
!     SBEND IS BASED ON SUBROUTINE BUILD OF THE SPARSE MATRIX       
!     PACKAGE SPARSPAK DEVELOPED BY ALAN GEORGE AND JOSEPH LUI      
!     OF THE UNIVERSITY OF WATERLOO, WATERLOO, ONTARIO.   
!       
! ... PARAMETERS  
!       
! ...... INPUT    
!       
!     N       THE ORDER OF THE LINEAR SYSTEM    
!       
!     NZ      THE LENGTH OF THE ARRAYS JA, IWORK, AND A.  
!       
! ...... INPUT/OUTPUT       
!       
!     IA      INTEGER ARRAY OF LENGTH N+1.  THE FIRST N ENTRIES     
!             POINT TO THE BEGINNING OF THE LINKED LIST FOR EACH    
!             ROW.  IA(N+1)-1 IS THE TOP OF THE LINKED LISTS
!             CONTAINED IN JA, IWORK, AND A.  ON OUTPUT IA WILL     
!             POINT TO THE FIRST ENTRY OF EACH ROW IN THE FINAL     
!             DATA STRUCTURE. 
!       
!     JA      INTEGER ARRAY OF LENGTH NZ.  ON INPUT JA STORES THE   
!             COLUMN NUMBERS OF THE NONZERO ENTRIES AS INDICATED    
!             BY THE LINKED LISTS.  ON OUTPUT JA STORES THE 
!             COLUMN NUMBERS IN ROW ORDERED FORM. 
!       
!     A       D.P. ARRAY OF LENGTH NZ.  ON INPUT A STORES THE       
!             VALUE OF THE NOZERO ENTRIES AS INDICATED BY THE       
!             LINKED LISTS.  ON OUTPUT A STORES THE VALUES IN       
!             ROW ORDERED FORM.       
!       
!     IWORK    INTEGER ARRAY OF LENGTH NZ.  ON INPUT IWORK STORES THE 
!             THE LINKS OF THE LINKED LISTS.  ON OUTPUT IT IS       
!             DESTROYED.    
!       
!***********************************************************************
!       
      INTEGER N,NZ,IA(1),JA(NZ),IWORK(NZ)       
      DOUBLE PRECISION A(NZ)
!       
      INTEGER MAXTOP,NEXT,TOP,IDEG,NULINK,JAJ,HLINK,OHLINK,L,I,LINK,&
     &   MHLINK   
      DOUBLE PRECISION VAL  
!       
!***********************************************************************
!       
! ... INITIALIZATION
!       
! ...... THE VARIABLES NEXT AND TOP RESPECTIVELY POINT TO THE       
!        NEXT AVAILABLE ENTRY FOR THE FINAL DATA STRUCTURE AND      
!        THE TOP OF THE REMAINDER OF THE LINKED LISTS.    
!       
      NEXT = 1    
      TOP = IA(N+1)+1       
      MAXTOP = NZ-IA(N+1)+1 
!       
!***********************************************************************
!       
! ... CONVERT EACH ROW INTO FINAL FORM
!       
      DO 90 I = 1,N 
         IDEG = 0 
         NULINK = IA(I)     
!       
! ... LOOP OVER EACH NODE IN THE LINKED LIST OF ROW I     
!       
   10    LINK = NULINK      
         IF (LINK.LE.0) GO TO 80      
         NULINK = IWORK(LINK) 
         JAJ = JA(LINK)     
         VAL = A(LINK)      
!       
! ... CHECK TO SEE IF A COLLISION BETWEEN THE LINKED LISTS
!     AND THE FINAL FORM HAS OCCURRED.
!       
         IF (NEXT.GE.TOP.AND.LINK.NE.TOP) GO TO 20
!       
! ... COLLISION HAS NOT OCCURRED.  FREE THE SPACE FOR THE TRIPLE    
!     (JA(LINK), A(LINK), IWORK(LINK))
!       
         JA(LINK) = 0       
         A(LINK) = 0.0D0    
         IWORK(LINK) = 0    
!       
! ... SPECIAL CASE TO MOVE  TOP  DOWN IF LINK .EQ. TOP    
!       
         IF (LINK.EQ.TOP) GO TO 60    
         GO TO 70 
!       
!***********************************************************************
!       
! ... COLLISION HAS OCCURRED.  CLEAR OFF SOME SPACE FOR THE CURRENT 
!     ENTRY BY MOVING THE TRIPLE ( JA(TOP),A(TOP),IWORK(TOP) )      
!     DOWNWARDS TO THE FREED TRIPLE ( JA(LINK),A(LINK),IWORK(LINK) ). 
!     THEN ADJUST THE LINK FIELDS.    
!       
! ...... PATCH UP THE LINKED LIST FOR THE CURRENT ROW I.  THEN      
!        TRAVERSE THE LINKED LIST CONTAINING TOP UNTIL THE POINTER  
!        POINTER BACK TO IA IS FOUND. 
!       
   20    IA(I) = LINK       
         HLINK = TOP
!       
   30    HLINK = IWORK(HLINK) 
         IF (HLINK.GT.0) GO TO 30     
!       
! ...... NOW FOLLOW THE LINKED LIST BACK TO TOP KEEPING TRACK       
!        OF THE OLD LINK.   
!       
! ......... SPECIAL CASE IF IA(-HLINK) = TOP    
!       
         MHLINK = -HLINK    
         IF (IA(MHLINK).NE.TOP) GO TO 40
!       
         IWORK(LINK) = IWORK(TOP)     
         JA(LINK) = JA(TOP) 
         A(LINK) = A(TOP)   
         IA(MHLINK) = LINK  
         IF (NULINK.EQ.TOP) NULINK = LINK       
         GO TO 60 
!       
! ......... USUAL CASE.     
!       
   40    HLINK = IA(MHLINK) 
   50    OHLINK = HLINK     
         HLINK = IWORK(OHLINK)
         IF (HLINK.NE.TOP) GO TO 50   
!       
         IWORK(LINK) = IWORK(TOP)     
         JA(LINK) = JA(TOP) 
         A(LINK) = A(TOP)   
         IF (OHLINK.NE.LINK) IWORK(OHLINK) = LINK 
         IF (NULINK.EQ.TOP) NULINK = LINK       
!       
! ... COLLAPSE TOP OF LINK LIST BY AS MUCH AS POSSIBLE    
!       
   60    TOP = TOP+1
         IF (TOP.GE.MAXTOP) GO TO 70  
         IF (IWORK(TOP).NE.0) GO TO 70
         GO TO 60 
!       
!***********************************************************************
!       
! ... PUT THE CURRENT TRIPLE INTO THE FINAL DATA STRUCTURE
!       
   70    JA(NEXT) = JAJ     
         A(NEXT) = VAL      
         NEXT = NEXT+1      
         IDEG = IDEG+1      
         GO TO 10 
!       
! ... FINAL STRUCTURE FOR ROW I IS COMPLETE.  LINKED LIST IS
!     DESTROYED AND WILL BE RECAPTURED AS NECESSARY BY THE
!     LOOP ON LABEL 60      
!       
   80    IA(I) = IDEG       
!       
   90 CONTINUE    
!       
!***********************************************************************
!       
! ... FINALIZE THE DATA STRUCTURE BY BUILDING THE FINAL VERSION OF  
!     IA. 
!       
      L = IA(1)+1 
      IA(1) = 1   
      DO 100 I = 1,N
         IDEG = IA(I+1)     
         IA(I+1) = L
         L = L+IDEG 
  100 CONTINUE    
!       
! ... FINAL IA, JA, A DATA STRUCTURE BUILT.     
!       
      RETURN      
      END











      SUBROUTINE SBINI (N,NZ,IA,JA,A,IWORK)     
!       
!***********************************************************************
!       
!     SBINI IS THE FIRST OF A SUITE OF THREE SUBROUTINES TO AID     
!     THE USER TO CONSTRUCT THE IA, JA, A DATA STRUCTURE USED       
!     IN ITPACK.  
!       
!     SBINI INITIALIZES THE ARRAYS IA, JA, IWORK, AND A.  THE OTHER 
!     SUBROUTINES IN THE SUITE ARE SBSIJ ( WHICH BUILDS A LINKED    
!     LIST REPRESENTATION OF THE MATRIX STRUCTURE ) AND SBEND ( WHICH 
!     RESTRUCTURE THE LINKED LIST FORM INTO THE FINAL FORM ).       
!       
! ... PARAMETERS  
!       
! ...... INPUT    
!       
!     N          THE ORDER OF THE LINEAR SYSTEM 
!       
!     NZ         THE MAXIMUM NUMBER OF NONZEROES ALLOWED IN THE     
!                LINEAR SYSTEM.       
!       
! ...... OUTPUT   
!       
!     IA         INTEGER ARRAY OF LENGTH N+1.  SBINI SETS THIS ARRAY
!                TO -I FOR I = 1 THRU N.  IA(N+1) IS SET TO NZ.     
!       
!     JA         INTEGER ARRAY OF LENGTH NZ.  INITIALIZED TO ZERO HERE. 
!       
!     A          D.P. ARRAY OF LENGTH NZ.  INITIALIZED TO ZERO HERE.
!       
!     IWORK       INTEGER ARRAY OF LENGTH NZ.  INITIALIZED TO ZERO HERE.
!       
!***********************************************************************
!       
      INTEGER N,NZ,IA(1),JA(NZ),IWORK(NZ),I     
      DOUBLE PRECISION A(NZ)
!       
!***********************************************************************
!       
      DO 10 I = 1,N 
         IA(I) = -I 
   10 CONTINUE    
      IA(N+1) = NZ
!       
      CALL IVFILL (NZ,JA,0) 
      CALL IVFILL (NZ,IWORK,0)
      CALL VFILL (NZ,A,0.D0)
!       
      RETURN      
      END












      SUBROUTINE SBSIJ (N,NZ,IA,JA,A,IWORK,II,JJ,VALL,MODE,LEVELL,NOUTT,&
     &   IERR)    
!       
!***********************************************************************
!       
!     SBSIJ IS THE SECOND OF A SUITE OF THREE SUBROUTINES TO AID IN 
!     THE CONSTRUCTION OF THE IA, JA, A DATA STRUCTURE USED IN      
!     ITPACK.     
!       
!     SBSIJ TAKES THE INDIVIDUAL ENTRIES OF THE SPARSE MATRIX AS    
!     GIVEN TO IT AT EACH CALL VIA  (I,J,VAL) AND INSERTS IT INTO   
!     A LINKED LIST REPRESENTATION OF THE SPARSE MATRIX.  
!       
!     EACH ROW OF THE SPARSE MATRIX IS ASSOCIATED WITH A CIRCULAR   
!     LINKED LIST BEGINNING AT IA(I).  THE LAST ENTERED ELEMENT IN  
!     EACH LIST POINTS BACK TO IA(I) WITH THE VALUE -I.  THE LINKS  
!     ARE STORED IN THE ARRAY IWORK, WHILE JA AND A STORE THE COLUMN
!     NUMBER AND VALUE IN PARALLEL TO IWORK.  THE LINKED LISTED ARE 
!     STORED BEGINNING AT ENTRY NZ AND WORKING BACKWARDS TOWARDS 1. 
!       
! ... PARAMETERS  
!       
! ...... INPUT    
!       
!     N       THE ORDER OF THE LINEAR SYSTEM    
!       
!     NZ      THE LENGTH OF THE ARRAYS  JA, A, AND IWORK  
!       
!     I, J    THE ROW AND COLUMN NUMBERS OF THE ENTRY OF THE SPARSE 
!             LINEAR SYSTEM TO BE ENTERED IN THE DATA STRUCTURE(=II,JJ) 
!       
!     VAL     THE NONZERO VALUE ASSOCIATED WITH (I,J)  (= VALL)     
!       
!     MODE    IF THE (I,J) ENTRY HAS ALREADY BEEN SET, MODE SPECIFIES 
!             THE WAY IN WHICH THE ENTRY IS TO BE TREATED.
!             IF   MODE .LT. 0  LET THE VALUE REMAIN AS IS
!                       .EQ. 0  RESET IT TO THE NEW VALUE 
!                       .GT. 0  ADD THE NEW VALUE TO THE OLD VALUE  
!       
!     NOUT  OUTPUT FILE NUMBER (= NOUTT)
!       
!     LEVEL   OUTPUT FILE SWITCH (= LEVELL)     
! ... INPUT/OUTPUT
!       
!     IA      INTEGER ARRAY OF LENGTH N+1.  THE FIRST N ENTRIES     
!             POINT TO THE BEGINNING OF THE LINKED LIST FOR EACH    
!             ROW.  IA(N+1) POINTS TO THE NEXT ENTRY AVAILABLE FOR  
!             STORING THE CURRENT ENTRY INTO THE LINKED LIST.       
!       
!     JA      INTEGER ARRAY OF LENGTH NZ.  JA STORES THE COLUMN     
!             NUMBERS OF THE NONZERO ENTRIES.   
!       
!     A       D.P. ARRAY OF LENGTH NZ.  A STORES THE VALUE OF THE   
!             NONZERO ENTRIES.
!       
!     IWORK   INTEGER ARRAY OF LENGTH NZ. IWORK STORES THE LINKS.   
!       
!     IER     ERROR FLAG.(= IERR)  POSSIBLE RETURNS ARE   
!             IER =    0   SUCCESSFUL COMPLETION
!                 =  700   ENTRY WAS ALREADY SET,  VALUE HANDLED    
!                          AS SPECIFIED BY MODE.
!                 =  701   IMPROPER VALUE OF EITHER I OR J INDEX    
!                 =  702   NO ROOM REMAINING, NZ TOO SMALL. 
!       
!***********************************************************************
!       
      INTEGER N,NZ,IA(1),JA(NZ),IWORK(NZ),II,JJ,MODE,LEVELL,NOUTT,IERR
      DOUBLE PRECISION A(NZ),VALL     
!       
      INTEGER LINK,NEXT,NPL1,I,J,LEVEL,NOUT,IER 
      DOUBLE PRECISION VAL,TEMP       
!       
!***********************************************************************
!       
! ... CHECK THE VALIDITY OF THE (I,J) ENTRY     
!       
      I = II      
      J = JJ      
      VAL = VALL  
      LEVEL = LEVELL
      NOUT = NOUTT
      IER = 0     
      IF (I.LE.0.OR.I.GT.N) IER = 701 
      IF (J.LE.0.OR.J.GT.N) IER = 701 
      IF (IER.EQ.0) GO TO 20
      IF (LEVEL.GE.0) WRITE (NOUT,10) IER,I,J   
   10 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    IN ITPACK ROUTINE SBSIJ   '/' ','    IER = ',I10/' ', &
     &   '    ( ',I10,' , ',I10,' )'/' ',       &
     &   '    IMPROPER VALUE FOR I OR J ')      
      GO TO 130   
!       
! ... TRAVERSE THE LINK LIST POINTED TO BY IA(I) UNTIL EITHER       
! ... THE J ENTRY OR THE END OF THE LIST HAS BEEN FOUND.  
!       
   20 NPL1 = N+1  
      LINK = IA(I)
!       
! ...... SPECIAL CASE FOR THE FIRST ENTRY IN THE ROW      
!       
      IF (LINK.GT.0) GO TO 30 
      NEXT = IA(NPL1)       
      IF (NEXT.LT.1) GO TO 110
!       
      IA(I) = NEXT
      JA(NEXT) = J
      A(NEXT) = VAL 
      IWORK(NEXT) = -I      
      IA(NPL1) = NEXT-1     
      GO TO 130   
!       
! ... FOLLOW THE LINK LIST UNTIL J OR THE END OF THE LIST IS FOUND  
!       
   30 IF (JA(LINK).EQ.J) GO TO 40     
      IF (IWORK(LINK).LE.0) GO TO 100 
      LINK = IWORK(LINK)    
      GO TO 30    
!       
!:      
! ... ENTRY (I,J) ALREADY HAS BEEN SET.  RESET VALUE DEPENDING ON MODE
!       
   40 IER = 700   
      IF (MODE.GE.0) GO TO 60 
      IF (LEVEL.GE.1) WRITE (NOUT,50) IER,I,J,A(LINK)     
   50 FORMAT ('0','*** W A R N I N G ************'/'0',   &
     &   '    IN ITPACK ROUTINE SBSIJ   '/' ','    IER = ',I10/' ', &
     &   '    ( ',I10,' , ',I10,' )'/' ',       &
     &   '    ENTRY ALREADY SET AND IS LEFT AS ',D15.8)   
      GO TO 130   
   60 IF (MODE.GE.1) GO TO 80 
      IF (LEVEL.GE.1) WRITE (NOUT,70) IER,I,J,A(LINK),VAL 
   70 FORMAT ('0','*** W A R N I N G ************'/'0',   &
     &   '    IN ITPACK ROUTINE SBSIJ   '/' ','    IER = ',I10/' ', &
     &   '    ( ',I10,' , ',I10,' )'/' ',       &
     &   '    ENTRY ALREADY SET - CURRENT VALUE OF',D15.8/' ',      &
     &   '                                RESET TO',D15.8)
      A(LINK) = VAL 
      GO TO 130   
   80 TEMP = A(LINK)+VAL    
      IF (LEVEL.GE.1) WRITE (NOUT,90) IER,I,J,A(LINK),TEMP
   90 FORMAT ('0','*** W A R N I N G ************'/'0',   &
     &   '    IN ITPACK ROUTINE SBSIJ   '/' ','    IER = ',I10/' ', &
     &   '    ( ',I10,' , ',I10,' )'/' ',       &
     &   '    ENTRY ALREADY SET - CURRENT VALUE OF',D15.8/' ',      &
     &   '                                RESET TO',D15.8)
      A(LINK) = TEMP
      GO TO 130   
!       
! ... ENTRY (I,J) HAS NOT BEEN SET.  ENTER IT INTO THE LINKED LIST  
!       
  100 NEXT = IA(NPL1)       
      IF (NEXT.LT.1) GO TO 110
!       
      IWORK(LINK) = NEXT    
      JA(NEXT) = J
      A(NEXT) = VAL 
      IWORK(NEXT) = -I      
      IA(NPL1) = NEXT-1     
      GO TO 130   
!       
!***********************************************************************
!       
! ... ERROR TRAP FOR NO ROOM REMAINING
!       
  110 IER = 702   
      IF (LEVEL.GE.0) WRITE (NOUT,120) IER      
  120 FORMAT ('0','*** F A T A L     E R R O R ************'/'0',   &
     &   '    IN ITPACK ROUTINE SBSIJ   '/' ','    IER = ',I10/' ', &
     &   '    NZ TOO SMALL - NO ROOM FOR NEW ENTRY')      
!       
  130 CONTINUE    
      IERR = IER  
      RETURN      
      END










      SUBROUTINE SCAL (NN,IA,JA,A,RHS,U,D,LEVEL,NOUT,IER) 
!       
! ... ORIGINAL MATRIX IS SCALED TO A UNIT DIAGONAL MATRIX.  RHS     
! ... AND U ARE SCALED ACCORDINGLY.  THE MATRIX IS THEN SPLIT AND   
! ... IA, JA, AND A RESHUFFLED.       
!       
! ... PARAMETER LIST:       
!       
!          N      DIMENSION OF MATRIX (= NN)    
!          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
!          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
!          RHS    RIGHT HAND SIDE OF MATRIX PROBLEM       
!          U      LATEST ESTIMATE OF SOLUTION   
!          D      OUTPUT VECTOR CONTAINING THE SQUARE ROOTS 
!                    OF THE DIAGONAL ENTRIES    
!          LEVEL  PRINTING SWITCH FOR ERROR CONDITION     
!          NOUT OUTPUT TAPE NUMBER    
!          IER    ERROR FLAG: ON RETURN NONZERO VALUES MEAN 
!                    401 : THE ITH DIAGONAL ELEMENT IS .LE. 0.      
!                    402 : NO DIAGONAL ELEMENT IN ROW I   
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!     
      INTEGER IA(NN+1),JA(3*NN),NN,LEVEL,NOUT,IER
!     INTEGER IA(1),JA(1),NN,LEVEL,NOUT,IER  //old version
      DOUBLE PRECISION A(3*NN),RHS(NN),U(NN),D(NN)
!     DOUBLE PRECISION A(1),RHS(NN),U(NN),D(NN)

! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER I,IBGN,IEND,II,IM1,J,JADD,JAJJ,JJ,JJPI,N,NP1
      DOUBLE PRECISION DI   
!       
! ... EXTRACT SQUARE ROOT OF THE DIAGONAL OUT OF A AND SCALE U AND RHS
!       
      N = NN      
      IER = 0     
      DO 80 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IF (IBGN.GT.IEND) GO TO 50   
         DO 40 JJ = IBGN,IEND 
            IF (JA(JJ).NE.II) GO TO 40
            DI = A(JJ)      
            IF (DI.GT.0.D0) GO TO 70  
            IF (DI.EQ.0.D0) GO TO 20  
            IER = 401       
            IF (LEVEL.GE.0) WRITE (NOUT,10) II  
   10       FORMAT ('0','*** F A T A L     E R R O R ************'/'0', &
     &         '    IN ITPACK ROUTINE SCAL    '/' ',      &
     &         '    DIAGONAL ENTRY IN ROW ',I10,' NEGATIVE')
            RETURN
   20       IER = 401       
            IF (LEVEL.GE.0) WRITE (NOUT,30)     
   30       FORMAT ('0','*** F A T A L     E R R O R ************'/'0', &
     &         '    IN ITPACK ROUTINE SCAL    '/' ',      &
     &         '    DIAGONAL ENTRY IN ROW ',I10,' IS ZERO') 
            RETURN
   40    CONTINUE 
   50    IER = 402
         IF (LEVEL.GE.0) WRITE (NOUT,60) II     
   60    FORMAT ('0','*** F A T A L     E R R O R ************'/'0',&
     &      '    IN ITPACK ROUTINE SCAL    '/' ', &
     &      '    NO DIAGONAL ENTRY IN ROW',I10) 
         RETURN   
!       
   70    CONTINUE 
         DI = DSQRT(DABS(DI)) 
         RHS(II) = RHS(II)/DI 
         U(II) = U(II)*DI   
         D(II) = DI 
   80 CONTINUE    
!       
! ... SHIFT MATRIX TO ELIMINATE DIAGONAL ENTRIES
!       
      IF (N.EQ.1) GO TO 110 
      NP1 = N+1   
      DO 100 I = 1,N
         IM1 = I-1
         II = NP1-I 
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         JADD = IBGN+IEND   
         DO 90 J = IBGN,IEND
            JJ = JADD-J     
            JJPI = JJ+IM1   
            IF (JA(JJ).EQ.II) IM1 = I 
            A(JJPI) = A(JJ) 
            JA(JJPI) = JA(JJ) 
   90    CONTINUE 
         IA(II+1) = IA(II+1)+I-1      
  100 CONTINUE    
  110 IA(1) = IA(1)+N       
!       
! ... SCALE SHIFTED MATRIX AND STORE D ARRAY IN FIRST N ENTRIES OF A
!       
      DO 140 II = 1,N       
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         DI = D(II) 
         IF (IBGN.GT.IEND) GO TO 130  
         DO 120 JJ = IBGN,IEND
            JAJJ = JA(JJ)   
            A(JJ) = A(JJ)/(DI*D(JAJJ))
  120    CONTINUE 
  130    CONTINUE 
         A(II) = DI 
  140 CONTINUE    
!       
      RETURN      
      END








      SUBROUTINE SUM3 (N,C1,X1,C2,X2,C3,X3)     
!       
! ... COMPUTES X3 = C1*X1 + C2*X2 + C3*X3       
!       
! ... PARAMETER LIST:       
!       
!          N        INTEGER LENGTH OF VECTORS X1, X2, X3  
!          C1,C2,C3 D.P. CONSTANTS    
!          X1,X2,X3 D.P. VECTORS SUCH THAT      
!                   X3(I) = C1*X1(I) + C2*X2(I) + C3*X3(I)
!                   X3(I) = C1*X1(I) + C2*X2(I)  IF C3 = 0. 
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER I   
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
      DOUBLE PRECISION X1(N),X2(N),X3(N),C1,C2,C3 
!       
      IF (N.LE.0) RETURN    
      IF (DABS(C3).EQ.0.D0) GO TO 20  
!       
      DO 10 I = 1,N 
         X3(I) = C1*X1(I)+C2*X2(I)+C3*X3(I)     
   10 CONTINUE    
      RETURN      
!       
! ... COMPUTE X3 = C1*X1 + C2*X2      
!       
   20 DO 30 I = 1,N 
         X3(I) = C1*X1(I)+C2*X2(I)    
   30 CONTINUE    
!       
      RETURN      
      END 
      DOUBLE PRECISION FUNCTION TAU (II)
!       
! ... THIS SUBROUTINE SETS TAU(II) FOR THE SOR METHOD.    
!       
! ... PARAMETER LIST:       
!       
!          II     NUMBER OF TIMES PARAMETERS HAVE BEEN CHANGED      
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER II  
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      DOUBLE PRECISION T(8) 
!       
      DATA T(1),T(2),T(3),T(4),T(5),T(6),T(7),T(8) / 1.5D0,1.8D0,1.85D0,&
     &   1.9D0,1.94D0,1.96D0,1.975D0,1.985D0 /  
!       
      TAU = 1.992D0 
      IF (II.LE.8) TAU = T(II)
!       
      RETURN      
      END 

      FUNCTION TIMER (TIMDMY) 
!
! ... TIMER IS A ROUTINE TO RETURN THE EXECUTION TIME IN
! ... SECONDS.
!
! ... PARAMETERS -- 
!
!          TIMDMY   DUMMY ARGUMENT
!
!
! *********************************************
! **                                         **
! **   THIS ROUTINE IS NOT PORTABLE.         **
! **                                         **
! *********************************************
!
      REAL TIMDMY
!
! ... CRAY Y-MP.
!
!     TIMER = SECOND ()
!
! ... UNIX ETIME FACILITY.
!
!      EXTERNAL ETIME
!      DIMENSION TARRAY(2)
!      REAL ETIME, TIMER
!      TOTAL = ETIME (TARRAY)
!      TIMER = TOTAL
!
! ... IBM RISC SYSTEM/6000.
!
!     TIMER = FLOAT(MCLOCK())/100.0
!

!     ZYL
      TIMER = 0
      RETURN
      END 







      LOGICAL FUNCTION TSTCHG (IBMTH) 
!       
!     THIS FUNCTION PERFORMS A TEST TO DETERMINE IF PARAMETERS      
!     SHOULD BE CHANGED FOR SEMI-ITERATION ACCELERATED METHODS.     
!       
! ... PARAMETER LIST:       
!       
!          IBMTH  INDICATOR OF BASIC METHOD BEING ACCELERATED BY SI 
!                      IBMTH = 1,   JACOBI      
!                            = 2,   REDUCED SYSTEM
!                            = 3,   SSOR
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IBMTH 
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER IP  
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCKS IN MAIN SUBROUTINE  
!       
      IP = IN-IS  
      IF (IBMTH.EQ.2) IP = 2*IP       
!       
      IF (IN.EQ.0) GO TO 10 
      IF (IP.LT.3) GO TO 20 
!       
      QA = DSQRT(DABS(DELNNM/DELSNM)) 
      QT = 2.D0*DSQRT(DABS(RRR**IP))/(1.D0+RRR**IP)       
      IF ((QA.GE.1.D0).OR.(QA.LT.QT**FF)) GO TO 20
!       
! ... TEST PASSES -- CHANGE PARAMETERS
!       
   10 TSTCHG = .TRUE.       
      RETURN      
!       
! ... TEST FAILS -- DO NOT CHANGE PARAMETERS    
!       
   20 TSTCHG = .FALSE.      
      RETURN      
!       
      END











      SUBROUTINE UNSCAL (N,IA,JA,A,RHS,U,D)     
!       
! ... THIS SUBROUTINE REVERSES THE PROCESS OF SCAL.       
!       
! ... PARAMETER LIST:       
!       
!          N      DIMENSION OF MATRIX 
!          IA,JA  INTEGER ARRAYS OF SPARSE MATRIX REPRESENTATION    
!          A      D.P. ARRAY OF SPARSE MATRIX REPRESENTATION
!          RHS    RIGHT HAND SIDE OF MATRIX PROBLEM       
!          U      LATEST ESTIMATE OF SOLUTION   
!          D      VECTOR CONTAINING THE SQUARE ROOTS      
!                    OF THE DIAGONAL ENTRIES    
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER IA(N+1),JA(3*N),N
      DOUBLE PRECISION A(3*N),RHS(N),U(N),D(N)
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER IBGN,IEND,II,INEW,IS,JAJJ,JJ,JJPI 
      DOUBLE PRECISION DI   
!       
! ... EXTRACT DIAGONAL FROM SCALED A AND UNSCALE U AND RHS
!       
      DO 10 II = 1,N
         DI = A(II) 
         U(II) = U(II)/DI   
         RHS(II) = RHS(II)*DI 
         D(II) = DI 
   10 CONTINUE    
!       
! ... UNSCALE A   
!       
      DO 30 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IF (IBGN.GT.IEND) GO TO 30   
         DI = D(II) 
         DO 20 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            A(JJ) = A(JJ)*DI*D(JAJJ)  
   20    CONTINUE 
   30 CONTINUE    
!       
! ... INSERT DIAGONAL BACK INTO A     
!       
      DO 60 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IS = N-II
         INEW = IBGN-IS-1   
         A(INEW) = D(II)**2 
         JA(INEW) = II      
         IF (IS.EQ.0.OR.IBGN.GT.IEND) GO TO 50  
         DO 40 JJ = IBGN,IEND 
            JJPI = JJ-IS    
            A(JJPI) = A(JJ) 
            JA(JJPI) = JA(JJ) 
   40    CONTINUE 
   50    CONTINUE 
         IA(II) = INEW      
   60 CONTINUE    
!       
      RETURN      
      END









      SUBROUTINE VEVMW (N,V,W)
!       
! ... VEVMW COMPUTES V = V - W
!       
! ... PARAMETER LIST:       
!       
!          N      INTEGER LENGTH OF VECTORS V AND W       
!          V      D.P. VECTOR 
!          W      D.P. VECTOR SUCH THAT   V(I) = V(I) - W(I)
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER N   
      DOUBLE PRECISION V(N),W(N)      
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER I,M,MP1       
!       
      IF (N.LE.0) RETURN    
      M = MOD(N,4)
!       
      IF (M.EQ.0) GO TO 20  
      DO 10 I = 1,M 
         V(I) = V(I)-W(I)   
   10 CONTINUE    
      IF (N.LT.4) RETURN    
!       
   20 MP1 = M+1   
      DO 30 I = MP1,N,4     
         V(I) = V(I)-W(I)   
         V(I+1) = V(I+1)-W(I+1)       
         V(I+2) = V(I+2)-W(I+2)       
         V(I+3) = V(I+3)-W(I+3)       
   30 CONTINUE    
      RETURN      
!       
      END











      SUBROUTINE VEVPW (N,V,W)
!       
! ... VPW COMPUTES    V = V + W       
!       
! ... PARAMETER LIST:       
!       
!          N      LENGTH OF VECTORS V AND W     
!          V      D.P. VECTOR 
!          W      D.P. VECTOR SUCH THAT   V(I) = V(I) + W(I)
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER N   
      DOUBLE PRECISION V(N),W(N)      
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER I,M,MP1       
!       
      IF (N.LE.0) RETURN    
!       
      M = MOD(N,4)
      IF (M.EQ.0) GO TO 20  
      DO 10 I = 1,M 
         V(I) = V(I)+W(I)   
   10 CONTINUE    
      IF (N.LT.4) RETURN    
!       
   20 MP1 = M+1   
      DO 30 I = MP1,N,4     
         V(I) = V(I)+W(I)   
         V(I+1) = V(I+1)+W(I+1)       
         V(I+2) = V(I+2)+W(I+2)       
         V(I+3) = V(I+3)+W(I+3)       
   30 CONTINUE    
!       
      RETURN      
      END










      SUBROUTINE VFILL (N,V,VAL)      
!       
!     FILLS A VECTOR, V, WITH A CONSTANT VALUE, VAL.      
!       
! ... PARAMETER LIST:       
!       
!          N      INTEGER LENGTH OF VECTOR V    
!          V      D.P. VECTOR 
!          VAL    D.P. CONSTANT THAT FILLS FIRST N LOCATIONS OF V   
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER N   
      DOUBLE PRECISION V(N),VAL       
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER I,M,MP1       
!       
      IF (N.LE.0) RETURN    
!       
!     CLEAN UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 10  
!       
      M = MOD(N,10) 
      IF (M.EQ.0) GO TO 20  
      DO 10 I = 1,M 
         V(I) = VAL 
   10 CONTINUE    
      IF (N.LT.10) RETURN   
!       
   20 MP1 = M+1   
      DO 30 I = MP1,N,10    
         V(I) = VAL 
         V(I+1) = VAL       
         V(I+2) = VAL       
         V(I+3) = VAL       
         V(I+4) = VAL       
         V(I+5) = VAL       
         V(I+6) = VAL       
         V(I+7) = VAL       
         V(I+8) = VAL       
         V(I+9) = VAL       
   30 CONTINUE    
!       
      RETURN      
      END












      SUBROUTINE VOUT (N,V,ISWT,NOUTT)
!       
!     THIS SUBROUTINE EFFECTS PRINTING OF RESIDUAL AND SOLUTION     
!     VECTORS - CALLED FROM PERROR5    
!       
! ... PARAMETER LIST:       
!       
!          V      VECTOR OF LENGTH N  
!          ISWT   LABELLING INFORMATION 
!          NOUT OUTPUT DEVICE NUMBER (= NOUTT)  
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER N,ISWT,NOUTT  
      DOUBLE PRECISION V(N) 
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER I,J,JM1,K,KUPPER,NOUT   
!       
      NOUT = NOUTT
!       
!        IF (N .LE. 0) RETURN 
!       
      KUPPER = MIN0(N,8)    
      IF (ISWT.EQ.1) WRITE (NOUT,10)  
   10 FORMAT (//5X,'RESIDUAL VECTOR') 
      IF (ISWT.EQ.2) WRITE (NOUT,20)  
   20 FORMAT (//5X,'SOLUTION VECTOR') 
      WRITE (NOUT,30) (I,I=1,KUPPER)  
   30 FORMAT (10X,8I15)     
      WRITE (NOUT,40)       
   40 FORMAT (10X,120('-')/)
!       
      DO 60 J = 1,N,8       
         KUPPER = MIN0(J+7,N) 
         JM1 = J-1
         WRITE (NOUT,50) JM1,(V(K),K=J,KUPPER)  
   50    FORMAT (4X,I5,'+  ',8D15.5)  
   60 CONTINUE    
!       
      RETURN      
      END










      SUBROUTINE WEVMW (N,V,W)
!       
! ... WEVMW COMPUTES W = V - W
!       
! ... PARAMETER LIST:       
!       
!          N      INTEGER LENGTH OF VECTORS V AND W       
!          V      D.P. VECTOR 
!          W      D.P. VECTOR SUCH THAT   W(I) = V(I) - W(I)
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      INTEGER N   
      DOUBLE PRECISION V(N),W(N)      
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      INTEGER I,M,MP1       
!       
      IF (N.LE.0) RETURN    
      M = MOD(N,4)
      IF (M.EQ.0) GO TO 20  
      DO 10 I = 1,M 
         W(I) = V(I)-W(I)   
   10 CONTINUE    
      IF (N.LT.4) RETURN    
!       
   20 MP1 = M+1   
      DO 30 I = MP1,N,4     
         W(I) = V(I)-W(I)   
         W(I+1) = V(I+1)-W(I+1)       
         W(I+2) = V(I+2)-W(I+2)       
         W(I+3) = V(I+3)-W(I+3)       
   30 CONTINUE    
!       
      RETURN      
      END









      SUBROUTINE ZBRENT (N,TRI,EPS,NSIG,AA,BB,MAXFNN,IER) 
!       
!   MODIFIED IMSL ROUTINE NAME   - ZBRENT       
!       
!-----------------------------------------------------------------------
!       
!   COMPUTER            - CDC/SINGLE  
!       
!   LATEST REVISION     - JANUARY 1, 1978       
!       
!   PURPOSE             - ZERO OF A FUNCTION WHICH CHANGES SIGN IN A
!                           GIVEN INTERVAL (BRENT ALGORITHM)
!       
!   USAGE               - CALL ZBRENT (F,EPS,NSIG,A,B,MAXFN,IER)    
!       
!   ARGUMENTS    TRI    - A TRIDIAGONAL MATRIX OF ORDER N 
!                EPS    - FIRST CONVERGENCE CRITERION (INPUT).  A ROOT, 
!                           B, IS ACCEPTED IF DABS(F(B)) IS LESS THAN OR
!                           EQUAL TO EPS.  EPS MAY BE SET TO ZERO.  
!                NSIG   - SECOND CONVERGENCE CRITERION (INPUT).  A ROOT,
!                           B, IS ACCEPTED IF THE CURRENT APPROXIMATION 
!                           AGREES WITH THE TRUE SOLUTION TO NSIG   
!                           SIGNIFICANT DIGITS. 
!                A,B    - ON INPUT, THE USER MUST SUPPLY TWO POINTS, A
!                           AND B, SUCH THAT F(A) AND F(B) ARE OPPOSITE 
!                           IN SIGN. (= AA, BB) 
!                           ON OUTPUT, BOTH A AND B ARE ALTERED.  B 
!                           WILL CONTAIN THE BEST APPROXIMATION TO THE
!                           ROOT OF F. SEE REMARK 1.      
!                MAXFN  - ON INPUT, MAXFN SHOULD CONTAIN AN UPPER BOUND 
!                           ON THE NUMBER OF FUNCTION EVALUATIONS   
!                           REQUIRED FOR CONVERGENCE.  ON OUTPUT, MAXFN 
!                           WILL CONTAIN THE ACTUAL NUMBER OF FUNCTION
!                           EVALUATIONS USED. (= MAXFNN)  
!                IER    - ERROR PARAMETER. (OUTPUT)       
!                         TERMINAL ERROR
!                           IER = 501 INDICATES THE ALGORITHM FAILED TO 
!                             CONVERGE IN MAXFN EVALUATIONS.
!                           IER = 502 INDICATES F(A) AND F(B) HAVE THE
!                             SAME SIGN.
!       
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32 
!                       - SINGLE/H36,H48,H60    
!       
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND       
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL  
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!       
!   REMARKS  1.  LET F(X) BE THE CHARACTERISTIC FUNCTION OF THE MATRIX
!                TRI EVALUATED AT X. FUNCTION DETERM EVALUATES F(X).
!                ON EXIT FROM ZBRENT, WHEN IER=0, A AND B SATISFY THE 
!                FOLLOWING, 
!                F(A)*F(B) .LE.0,     
!                DABS(F(B)) .LE. DABS(F(A)), AND
!                EITHER DABS(F(B)) .LE. EPS OR  
!                DABS(A-B) .LE. MAX(DABS(B),0.1)*10.0**(-NSIG).     
!                THE PRESENCE OF 0.1 IN THIS ERROR CRITERION CAUSES 
!                LEADING ZEROES TO THE RIGHT OF THE DECIMAL POINT TO BE 
!                COUNTED AS SIGNIFICANT DIGITS. SCALING MAY BE REQUIRED 
!                IN ORDER TO ACCURATELY DETERMINE A ZERO OF SMALL   
!                MAGNITUDE. 
!            2.  ZBRENT IS GUARANTEED TO REACH CONVERGENCE WITHIN   
!                K = (DLOG((B-A)/D)+1.0)**2 FUNCTION EVALUATIONS WHERE
!                  D=MIN(OVER X IN (A,B) OF     
!                    MAX(DABS(X),0.1)*10.0**(-NSIG)).     
!                THIS IS AN UPPER BOUND ON THE NUMBER OF EVALUATIONS. 
!                RARELY DOES THE ACTUAL NUMBER OF EVALUATIONS USED BY 
!                ZBRENT EXCEED DSQRT(K). D CAN BE COMPUTED AS FOLLOWS,
!                  P = DBLE(AMIN1(DABS(A),DABS(B)))       
!                  P = DMAX1(0.1,P)   
!                  IF ((A-0.1)*(B-0.1).LT.0.0) P = 0.1    
!                  D = P*10.0**(-NSIG)
!       
!   COPYRIGHT           - 1977 BY IMSL, INC. ALL RIGHTS RESERVED.   
!       
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.    
!       
!-----------------------------------------------------------------------
!       
! *** BEGIN: ITPACK COMMON  
!       
      INTEGER IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      COMMON /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      LOGICAL ADAPT,BETADT,CASEII,HALT,PARTAD   
      COMMON /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
      DOUBLE PRECISION BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA,&
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
      COMMON /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, &
     &   QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA      
!       
! *** END  : ITPACK COMMON  
!       
!     DESCRIPTION OF VARIABLES IN COMMON BLOCK IN MAIN SUBROUTINE   
!       
!                                  SPECIFICATIONS FOR ARGUMENTS     
!       
      INTEGER NSIG,MAXFNN,IER 
      DOUBLE PRECISION TRI(2,1),EPS,AA,BB       
!       
!                                  SPECIFICATIONS FOR LOCAL VARIABLES 
!       
      INTEGER IC,MAXFN      
      DOUBLE PRECISION ZERO,HALF,ONE,THREE,TEN,A,B,T,FA,FB,C,FC,D,E,TOL,&
     &   RM,S,P,Q,R,RONE,TEMP,DETERM  
      DATA ZERO / 0.D0 / ,HALF / 5.D-1 / ,ONE / 1.D0 / ,THREE / 3.D0 / ,&
     &   TEN / 10.D0 /      
!       
!                                  FIRST EXECUTABLE STATEMENT       
!       
      A = AA      
      B = BB      
      MAXFN = MAXFNN
      IER = 0     
      T = TEN**(-NSIG)      
      IC = 2      
      FA = DETERM(N,TRI,A)  
      FB = DETERM(N,TRI,B)  
      S = B       
!       
!                                  TEST FOR SAME SIGN     
!       
      IF (FA*FB.GT.ZERO) GO TO 110    
   10 C = A       
      FC = FA     
      D = B-C     
      E = D       
   20 IF (DABS(FC).GE.DABS(FB)) GO TO 30
      A = B       
      B = C       
      C = A       
      FA = FB     
      FB = FC     
      FC = FA     
   30 CONTINUE    
      TOL = T*DMAX1(DABS(B),0.1D0)    
      RM = (C-B)*HALF       
!       
!                                  TEST FOR FIRST CONVERGENCE CRITERIA
!       
      IF (DABS(FB).LE.EPS) GO TO 80   
!       
!                                  TEST FOR SECOND CONVERGENCE CRITERIA 
!       
      IF (DABS(C-B).LE.TOL) GO TO 80  
!       
!                                  CHECK EVALUATION COUNTER 
!       
      IF (IC.GE.MAXFN) GO TO 90       
!       
!                                  IS BISECTION FORCED    
!       
      IF (DABS(E).LT.TOL) GO TO 60    
      IF (DABS(FA).LE.DABS(FB)) GO TO 60
      S = FB/FA   
      IF (A.NE.C) GO TO 40  
!       
!                                  LINEAR INTERPOLATION   
!       
      P = (C-B)*S 
      Q = ONE-S   
      GO TO 50    
!       
!                                  INVERSE QUADRATIC INTERPOLATION  
!       
   40 Q = FA/FC   
      R = FB/FC   
      RONE = R-ONE
      P = S*((C-B)*Q*(Q-R)-(B-A)*RONE)
      Q = (Q-ONE)*RONE*(S-ONE)
   50 IF (P.GT.ZERO) Q = -Q 
      IF (P.LT.ZERO) P = -P 
      S = E       
      E = D       
!       
!                                  IF DABS(P/Q).GE.75*DABS(C-B) THEN
!                                     FORCE BISECTION     
!       
      IF (P+P.GE.THREE*RM*Q) GO TO 60 
!       
!                                  IF DABS(P/Q).GE..5*DABS(S) THEN FORCE
!                                     BISECTION. S = THE VALUE OF P/Q 
!                                     ON THE STEP BEFORE THE LAST ONE 
!       
      IF (P+P.GE.DABS(S*Q)) GO TO 60  
      D = P/Q     
      GO TO 70    
!       
!                                  BISECTION    
!       
   60 E = RM      
      D = E       
!       
!                                  INCREMENT B  
!       
   70 A = B       
      FA = FB     
      TEMP = D    
      IF (DABS(TEMP).LE.HALF*TOL) TEMP = DSIGN(HALF*TOL,RM) 
      B = B+TEMP  
      S = B       
      FB = DETERM(N,TRI,S)  
      IC = IC+1   
      IF (FB*FC.LE.ZERO) GO TO 20     
      GO TO 10    
!       
!                                  CONVERGENCE OF B       
!       
   80 A = C       
      MAXFN = IC  
      GO TO 130   
!       
!                                  MAXFN EVALUATIONS      
!       
   90 IER = 501   
      A = C       
      MAXFN = IC  
      IF (LEVEL.GE.1) WRITE (NOUT,100) MAXFN    
  100 FORMAT ('0','*** W A R N I N G ************'/'0',   &
     &   '    IN ITPACK ROUTINE ZBRENT  '/' ',  &
     &   '    ALGORITHM FAILED TO CONVERGE   '/' ','    IN',I6,     &
     &   ' ITERATIONS ')    
      GO TO 130   
!       
!                                  TERMINAL ERROR - F(A) AND F(B) HAVE
!                                  THE SAME SIGN
!       
  110 IER = 502   
      MAXFN = IC  
      IF (LEVEL.GE.1) WRITE (NOUT,120)
  120 FORMAT ('0','*** W A R N I N G ************'/'0',   &
     &   '    IN ITPACK ROUTINE ZBRENT  '/' ',  &
     &   '    F(A) AND F(B) HAVE SAME SIGN   ') 
  130 CONTINUE    
      AA = A      
      BB = B      
      MAXFNN = MAXFN
      RETURN      
      END 
