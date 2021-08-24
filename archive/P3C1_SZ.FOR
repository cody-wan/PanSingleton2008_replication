C  SZ START MAIN ========================================================


C
C    THIS PROGRAM ESTIMATES
C    ONE FACTOR PAN/SINGLETON
C    MODEL VIA MLE
C
C    BULGARIA


C  SZ Comment out ISMLF77
C     USE IMSLF77
      IMPLICIT REAL*8 (A-H,P-Z)
      EXTERNAL FCN

C  SZ ADJUST DATA, PFILE AS NECESSARY
C  SZ 260 is from TPFILE, 260 COLUMNS FOR 1W to 260W MATURITIES
C  SZ NOTE DISCOUNT FUNCTION IS ONLY 2003m1 TO 2010M2
      DIMENSION PARM(10)
      DIMENSION DATA(85,15,3)
      DIMENSION PFILE(85,260)
      DIMENSION PIN(7)
      DIMENSION XGUESS(7)
      DIMENSION XLB(7),XUB(7),SOLV(7)
      DIMENSION ICAR(10)
      DIMENSION XSCALE(7),IPARAM(7),RPARAM(7)

      COMMON ICAR,PARM,DATA,PFILE

C  ADJUST FILE NAMES AS NECESSARY
      OPEN(UNIT=69,FILE="CDS1Y.TXT",STATUS='UNKNOWN')
      OPEN(UNIT=11,FILE="CDS3Y.TXT",STATUS='UNKNOWN')
      OPEN(UNIT=12,FILE="CDS5Y.TXT",STATUS='UNKNOWN')
      OPEN(UNIT=13,FILE='TPFILE.',STATUS='UNKNOWN')

C  SZ ADJUST OUTPUT NAME
      OPEN(UNIT=20,FILE='OUT1.',STATUS='UNKNOWN')

C  SZ Set default printing/stopping for errors
      IERSVR = 0
      IPACT  = 0
      ISACT  = 0
      CALL ERSET(IERSVR,IPACT,ISACT)

C
C  READ IN DATA
C  SZ ADJUST 86 FOR DATA and 16 for COUNTRIES AS NECESSARY
C  SZ I indexes rows (dates), J indexes columns (countries),
C  SZ and 1,2,3 for CDS 1, 3 and 5 year data 
      DO 10 I=1,85
      READ(69,*)(DATA(I,J,1),J=1,15)
      READ(11,*)(DATA(I,J,2),J=1,15)
      READ(12,*)(DATA(I,J,3),J=1,15)
      READ(13,*)(PFILE(I,J),J=1,260)
10    CONTINUE

      DO 15 I=1,85
      DO 15 J=1,15
      DATA(I,J,1)=DATA(I,J,1)/10000.
      DATA(I,J,2)=DATA(I,J,2)/10000.
      DATA(I,J,3)=DATA(I,J,3)/10000.
15    CONTINUE

C
C  PARAMETERS
C   1: ALPHA
C   2: BETA
C   3: SIGMA
C   4: ALPHA HAT
C   5: BETA HAT
C   6: SIG ERR1
C   7: SIG ERR2
      N=7

      xguess(1) =  3.520055
      xguess(2) = -0.8266322
      xguess(3)  =  .7495782
      XGUESS(4)=  -3.3802954
      XGUESS(5)=    .6875973
      XGUESS(6)=    .005881
      XGUESS(7)=    .003547

      IBTYPE=0

      XLB(1)=  -7.
      XLB(2)=  -2.
      XLB(3)=  .0001
      XLB(4)=  -7.
      XLB(5)=  .0001
      XLB(6)=  .00001
      XLB(7)=  .00001

      XUB(1)=  7.
      XUB(2)=   .5
      XUB(3)=  5.
      XUB(4)=   10.
      XUB(5)=   10.
      XUB(6)= .0500
      XUB(7)= .0500

C  SZ scaling 1 for 7 params
      DO 30 IJK=1,7
      XSCALE(IJK)=1.
30    CONTINUE

C  SZ call FCN 
      FTOL=.00000000000000001
      MAXFCN=5000
      IPARAM(1)=0
      FSCALE = 1000.

c      CALL DBCPOL(FCN,N,XGUESS,IBTYPE,XLB,XUB,FTOL,
c     * MAXFCN,SOLV,FVALUE)

c      CALL DBCONF(FCN,N,XGUESS,
c     *            IBTYPE,XLB,XUB,XSCALE,FSCALE,
c     *            IPARAM,RPARAM,SOLV,FVALUE)

      CALL FCN(N,XGUESS,FVALUE)



100   CONTINUE

1000  FORMAT(3X,60F12.7)
      STOP
      END

C  SZ END MAIN ========================================================





C  SZ START FCN ========================================================


      SUBROUTINE FCN(N,PIN,FVAL)
      IMPLICIT REAL*8 (A-H,P-Z)


C  SZ ADJUST DATA, PFILE DIMENSIONS
      DIMENSION PARM(10)
      DIMENSION DATA(85,15,3)
      DIMENSION PFILE(85,260)
      DIMENSION ICAR(10)
      DIMENSION PIN(7)
      DIMENSION CDS1T(200),CDS3T(200),CDS5T(200)
      DIMENSION CDS1TP(200),CDS3TP(200),CDS5TP(200)
      DIMENSION X(200)

      COMMON ICAR,PARM,DATA,PFILE

      PARM(1)  = PIN(1)
      PARM(2)  = PIN(2)
      PARM(3)  = PIN(3)
      PARM(4)  = PIN(4)
      PARM(5)  = PIN(5)
      PARM(6)  = PIN(6)
      PARM(7)  = PIN(7)

      ALP   = PARM(1)
      BET   = PARM(2)
      SIG   = PARM(3)
      ALPH  = PARM(4)
      BETH  = PARM(5)
      SIG1  = PARM(6)
      SIG2  = PARM(7)

      IF(SIG.LE.0.)SIG=.0001
      IF(SIG1.LE.0.)SIG1=.00001
      IF(SIG2.LE.0.)SIG2=.00001
      IF(BETH.LE.0.)BETH=.001

C
C  SWITCH FOR COUNTRY
C     INDICATOR
C
C     ICAR(1)=1
      ICAR(1)=1
      IN = ICAR(1)


C
C  LOOP THROUGH, SOLVE FOR
C  XLAM, AND THEN EVALUATE
C  FIT ACROSS ALL CDS CONTRACTS
C

C  SZ CHANGE 86 LOOP AS NECESSARY
      TOTAL = 0.
      TOTALL = 0.

      DO 100 I=1,85

C
C  SOLVE FOR
C   ALL CDS PRICES
C   GIVEN PARAMETER
C   VECTOR
C

      ICAR(2)=I

      CALL EVAL(X,CDS1T,CDS3T,CDS5T)


      CDS1 = DATA(I,IN,1)
      CDS3 = DATA(I,IN,2)
      CDS5 = DATA(I,IN,3)

C
C  SOLVE FOR XLAM
C

      CALL XINTP(CDS3,X,XLAM,CDS1T,CDS3T,
     *           CDS5T,CDS1X,CDS3X,CDS5X,DERIV)

      TOTALL = TOTALL + (CDS1-CDS1X)**2.
     *              + (CDS3-CDS3X)**2.
     *              + (CDS5-CDS5X)**2.


C
C  SOLVE FOR VALUES OF CDS
C   UNDER P MEASURE PARAMETERS
C

      CALL EVALP(X,CDS1TP,CDS3TP,CDS5TP)
      CALL XINTPP(CDS3,X,XLAM,CDS1TP,CDS3TP,
     *           CDS5TP,CDS1XP,CDS3XP,CDS5XP)


C  SZ OUTPUT file is written here, and here is the order
C  SZ XP is default risk
      write(20,1000)XLAM,CDS1,CDS1X,CDS1XP,
     *                   CDS3,CDS3X,CDS3XP,
     *                   CDS5,CDS5X,CDS5XP


       IF(I.EQ.1) GO TO 90

C
C  SOLVE FOR LOG LIKELIHOOD
C
      ER1 = CDS1-CDS1X
      ER2 = CDS5-CDS5X

      V = SIG*SIG/(2.*BETH)*(1.-EXP(-2.*BETH/12.))
      XM =LOG(XLAMM)*EXP(-BETH/12.) +ALPH/BETH *
     *                 (1.-EXP(-BETH/12.))

      PI=3.14159265


C      WRITE(*,1000)SIG1,SIG2,V,XLAM,XLAMM,DERIV

      TOTAL = TOTAL-.5*LOG(2.*PI*SIG1**2.)
     *             -ER1*ER1/(2.*SIG1**2.)
     *             -.5*LOG(2.*PI*SIG2**2.)
     *             -ER2*ER2/(2.*SIG2**2.)
     *             -.5*LOG(2.*PI*V)
     *             -(LOG(XLAM)-XM)**2./(2.*V)
     *             -LOG(DERIV)

90    XLAMM=XLAM

100   CONTINUE

C  SZ CHANGE this 86???
C      FVALL=TOTALL/(3.*86.)
      FVALL=TOTALL/(3.*85.)
      FVALL=FVALL**.5

      FVAL = -TOTAL

c      FVAL = FVALL

c      FVAL = -TOTAL+100000.*FVALL


C  SZ first write just outputs params
C  SZ second write outputs params in last row of output
C  SZ third is debuggin country column
      WRITE(*,1000)pin(1),pin(2),pin(3),FVAL,FVALL
      write(20,1000)(pin(j),j=1,7),fval,FVALL
      WRITE(*,1002)IN

1000  FORMAT(3X,100F16.7)
1001  FORMAT(3X,10F16.5)
1002  FORMAT(I4,F8.3)
      RETURN
      END

C  SZ END OF FCN ========================================================








C  SZ START EVAL ========================================================


      SUBROUTINE EVAL(X,CDS1T,CDS3T,CDS5T)
      IMPLICIT REAL*8 (A-H,P-Z)

C
C    THIS SUBROUTINE EVALUATES
C     THE MODEL VALUES OF CDS1,
C     CDS3, AND CDS5
C

C  SZ CHANGE DATA, PFILE
      DIMENSION PARM(10)
      DIMENSION DATA(85,15,3)
      DIMENSION PFILE(85,260)
      DIMENSION ICAR(10)
      DIMENSION X(200)
      DIMENSION F(200,61)
      DIMENSION G(200,61)
      DIMENSION CDS1T(200),CDS3T(200),CDS5T(200)

      COMMON ICAR,PARM,DATA,PFILE

      WR=.75

      IN=ICAR(1)
      II=ICAR(2)

      CALL PDE1(PARM,X,F)
      CALL PDE2(PARM,X,G)

      DO 10 I=1,200

C
C  VALUE 1 YR CDS
C

      SUM1=0.
      SUM2=0.
      DO 100 J=1,4
      SUM1=SUM1+PFILE(II,13*J)*F(I,1+J*3)
      SUM2=SUM2+PFILE(II,13*J)*G(I,1+J*3)
100   CONTINUE
      CDS1T(I) = WR*SUM2/SUM1

C
C  VALUE 3 YR CDS
C

      SUM1=0.
      SUM2=0.
      DO 110 J=1,12
      SUM1=SUM1+PFILE(II,13*J)*F(I,1+J*3)
      SUM2=SUM2+PFILE(II,13*J)*G(I,1+J*3)
110   CONTINUE
      CDS3T(I) = WR*SUM2/SUM1

C
C  VALUE 5 YR CDS
C

      SUM1=0.
      SUM2=0.
      DO 120 J=1,20
      SUM1=SUM1+PFILE(II,13*J)*F(I,1+J*3)
      SUM2=SUM2+PFILE(II,13*J)*G(I,1+J*3)
120   CONTINUE
      CDS5T(I) = WR*SUM2/SUM1

10    CONTINUE


1000  FORMAT(3X,100F16.7)

      RETURN
      END

C  SZ END EVAL ========================================================








C  SZ START PDE1 ========================================================


      SUBROUTINE PDE1(PARM,XX,F)
      IMPLICIT REAL*8 (A-H,P-Z)

C
C  THIS SUBROUTINE SOLVES THE
C  PDE FOR THE FIRST TYPE OF
C  INITIAL CONDITION AND RETURNS
C  THE MATRIX OF VALUES AS F
C

      DIMENSION X(200),XX(200)
      DIMENSION F(200,61)
      DIMENSION A(200),B(200,200)
      DIMENSION PARM(10)
      DIMENSION C1(200),C2(200),C3(200)

C
C  NEED TO SET LOWER LIMIT
C   AND STEP SIZES
C

      XLOWER = .0001
      XUPPER = .2500

      TSTEP=  1./12.
      XSTEP=  (LOG(XUPPER)-LOG(XLOWER))/200.

      ALP  = PARM(1)
      BET  = PARM(2)
      SIG  = PARM(3)
      ALPH = PARM(4)
      BETH = PARM(5)

      DO 10 I=1,200
      X(I)=LOG(XLOWER)+XSTEP*FLOAT(I-1)
      XX(I)=EXP(X(I))
C
C  SWITCH FOR INITIAL CONDITION
C

      F(I,1)=1.
C      F(I,1)=EXP(X(I))

10    CONTINUE

      DO 100 IT=2,61

      DO 110 I=1,200
      DO 110 J=1,200
      B(I,J)=0.
110   CONTINUE

      DO 120 I=1,200
      A(I) =   - F(I,IT-1)/TSTEP
120   CONTINUE

      B(1,1)= -(ALP-BET*X(1))/XSTEP-EXP(X(1))-1./TSTEP
      B(1,2)= (ALP-BET*X(1))/XSTEP
      B(200,199)= -(ALP-BET*X(100))/XSTEP
      B(200,200)= (ALP-BET*X(100))/XSTEP-EXP(X(100))-1./TSTEP

      DO 130 I=2,199
      B(I,I-1)=  SIG**2./(2.*XSTEP**2.)-(ALP-BET*X(I))/(2.*XSTEP)
      B(I,I)=   -SIG**2./XSTEP**2.-EXP(X(I))-1./TSTEP
      B(I,I+1)=  SIG**2./(2.*XSTEP**2.)+(ALP-BET*X(I))/(2.*XSTEP)
130   CONTINUE

      DO 140 I=1,199
      C1(I)=B(I,I+1)
140   CONTINUE

      DO 150 I=1,200
      C2(I)=B(I,I)
150   CONTINUE

      DO 160 I=2,200
      C3(I)=B(I,I-1)
160   CONTINUE

      N=200
      CALL DLSLTR(N,C3,C2,C1,A)

      DO 170 I=1,200
      F(I,IT)=A(I)
170   CONTINUE

100   CONTINUE


1000  FORMAT(3X,50F12.7)

      END

C  SZ END PDE1 ========================================================






C  SZ START PDE2 ========================================================

      SUBROUTINE PDE2(PARM,XX,G)
      IMPLICIT REAL*8 (A-H,P-Z)

C
C  THIS SUBROUTINE SOLVES THE
C  PDE FOR THE SECOND TYPE OF
C  INITIAL CONDITION AND RETURNS
C  THE MATRIX OF VALUES AS F
C

      DIMENSION X(200),XX(200)
      DIMENSION G(200,61)
      DIMENSION A(200),B(200,200)
      DIMENSION PARM(10)
      DIMENSION C1(200),C2(200),C3(200)

C
C  NEED TO SET LOWER LIMIT
C   AND STEP SIZES
C

      XLOWER = .0001
      XUPPER = .2500

      TSTEP=  1./12.
      XSTEP=  (LOG(XUPPER)-LOG(XLOWER))/200.

      ALP  = PARM(1)
      BET  = PARM(2)
      SIG  = PARM(3)
      ALPH = PARM(4)
      BETH = PARM(5)

      DO 10 I=1,200
      X(I)=LOG(XLOWER)+XSTEP*FLOAT(I-1)
      XX(I)=EXP(X(I))

C
C  SWITCH FOR INITIAL CONDITION
C

C      G(I,1)=1.
      G(I,1)=EXP(X(I))

10    CONTINUE

      DO 100 IT=2,61

      DO 110 I=1,200
      DO 110 J=1,200
      B(I,J)=0.
110   CONTINUE

      DO 120 I=1,200
      A(I) =   - G(I,IT-1)/TSTEP
120   CONTINUE

      B(1,1)= -(ALP-BET*X(1))/XSTEP-EXP(X(1))-1./TSTEP
      B(1,2)= (ALP-BET*X(1))/XSTEP
      B(200,199)= -(ALP-BET*X(100))/XSTEP
      B(200,200)= (ALP-BET*X(100))/XSTEP-EXP(X(100))-1./TSTEP

      DO 130 I=2,199
      B(I,I-1)=  SIG**2./(2.*XSTEP**2.)-(ALP-BET*X(I))/(2.*XSTEP)
      B(I,I)=   -SIG**2./XSTEP**2.-EXP(X(I))-1./TSTEP
      B(I,I+1)=  SIG**2./(2.*XSTEP**2.)+(ALP-BET*X(I))/(2.*XSTEP)
130   CONTINUE

      DO 140 I=1,199
      C1(I)=B(I,I+1)
140   CONTINUE

      DO 150 I=1,200
      C2(I)=B(I,I)
150   CONTINUE

      DO 160 I=2,200
      C3(I)=B(I,I-1)
160   CONTINUE

      N=200
      CALL DLSLTR(N,C3,C2,C1,A)

      DO 170 I=1,200
      G(I,IT)=A(I)
170   CONTINUE

100   CONTINUE


1000  FORMAT(3X,50F12.7)

      END

C  SZ END PDE2 =======================================================






C  SZ START XINTP ========================================================

      SUBROUTINE XINTP(CDS3,X,XLAM,
     *                 CDS1T,CDS3T,CDS5T,
     *                 CDS1X,CDS3X,CDS5X,DERIV)

      IMPLICIT REAL*8 (A-H,P-Z)

C
C    THIS SUBROUTINE INTERPOLATES
C    MATRIX OF CDS PRICES TO MATCH
C    MARKET DATA
C

C  SZ CHANGE DATA, PFILE
      DIMENSION PARM(10)
      DIMENSION DATA(85,15,3)
      DIMENSION PFILE(85,260)
      DIMENSION ICAR(10)
      DIMENSION X(200)
      DIMENSION F(200,61)
      DIMENSION G(200,61)
      DIMENSION CDS1T(200),CDS3T(200),CDS5T(200)

      COMMON ICAR,PARM,DATA,PFILE

      DO 10 I=1,199
      IF((CDS3.GE.CDS3T(I)).AND.(CDS3.LE.CDS3T(I+1)))THEN

      WEIGHT = (CDS3-CDS3T(I))/(CDS3T(I+1)-CDS3T(I))
      XLAM   = (1.-WEIGHT)*X(I)+WEIGHT*X(I+1)
      CDS1X  = (1.-WEIGHT)*CDS1T(I)+WEIGHT*CDS1T(I+1)
      CDS3X  = (1.-WEIGHT)*CDS3T(I)+WEIGHT*CDS3T(I+1)
      CDS5X  = (1.-WEIGHT)*CDS5T(I)+WEIGHT*CDS5T(I+1)

      DERIV =  (CDS3T(I+1)-CDS3T(I))/(LOG(X(I+1))
     *         -LOG(X(I)))


C      WRITE(*,1000)CDS3T(I+1),CDS3T(I),LOG(X(I+1)),LOG(X(I))

      GO TO 20
      ENDIF

10    CONTINUE
      XLAM=X(200)
      CDS1X=CDS1T(200)
      CDS3X=CDS3T(200)
      CDS5X=CDS5T(200)


1000  FORMAT(3X,20F16.7)

20    CONTINUE
      RETURN
      END

C  SZ END XINTP ========================================================





C  SZ START EVALP ========================================================

      SUBROUTINE EVALP(X,CDS1T,CDS3T,CDS5T)
      IMPLICIT REAL*8 (A-H,P-Z)

C
C    THIS SUBROUTINE EVALUATES
C     THE MODEL VALUES OF CDS1,
C     CDS3, AND CDS5
C
C     UNDER P MEASURE
C

C  SZ CHANGE DATA, PFILE

      DIMENSION PARM(10)
      DIMENSION DATA(85,15,3)
      DIMENSION PFILE(85,260)
      DIMENSION ICAR(10)
      DIMENSION X(200)
      DIMENSION F(200,61)
      DIMENSION G(200,61)
      DIMENSION CDS1T(200),CDS3T(200),CDS5T(200)

      COMMON ICAR,PARM,DATA,PFILE

      WR=.75

      IN=ICAR(1)
      II=ICAR(2)

      CALL PDE1P(PARM,X,F)
      CALL PDE2P(PARM,X,G)

      DO 10 I=1,200

C
C  VALUE 1 YR CDS
C

      SUM1=0.
      SUM2=0.
      DO 100 J=1,4
      SUM1=SUM1+PFILE(II,13*J)*F(I,1+J*3)
      SUM2=SUM2+PFILE(II,13*J)*G(I,1+J*3)
100   CONTINUE
      CDS1T(I) = WR*SUM2/SUM1

C
C  VALUE 3 YR CDS
C

      SUM1=0.
      SUM2=0.
      DO 110 J=1,12
      SUM1=SUM1+PFILE(II,13*J)*F(I,1+J*3)
      SUM2=SUM2+PFILE(II,13*J)*G(I,1+J*3)
110   CONTINUE
      CDS3T(I) = WR*SUM2/SUM1

C
C  VALUE 5 YR CDS
C

      SUM1=0.
      SUM2=0.
      DO 120 J=1,20
      SUM1=SUM1+PFILE(II,13*J)*F(I,1+J*3)
      SUM2=SUM2+PFILE(II,13*J)*G(I,1+J*3)
120   CONTINUE
      CDS5T(I) = WR*SUM2/SUM1

10    CONTINUE


1000  FORMAT(3X,100F16.7)

      RETURN
      END

C  SZ END EVALP ========================================================





C  SZ START PDE1P ========================================================


      SUBROUTINE PDE1P(PARM,XX,F)
      IMPLICIT REAL*8 (A-H,P-Z)

C
C  THIS SUBROUTINE SOLVES THE
C  PDE FOR THE FIRST TYPE OF
C  INITIAL CONDITION AND RETURNS
C  THE MATRIX OF VALUES AS F
C
C  USES P MEASURE
C

      DIMENSION X(200),XX(200)
      DIMENSION F(200,61)
      DIMENSION A(200),B(200,200)
      DIMENSION PARM(10)
      DIMENSION C1(200),C2(200),C3(200)

C
C  NEED TO SET LOWER LIMIT
C   AND STEP SIZES
C

      XLOWER = .0001
      XUPPER = .2500

      TSTEP=  1./12.
      XSTEP=  (LOG(XUPPER)-LOG(XLOWER))/200.

      ALP  = PARM(1)
      BET  = PARM(2)
      SIG  = PARM(3)
      ALPH = PARM(4)
      BETH = PARM(5)

      ALP=ALPH
      BET=BETH


      DO 10 I=1,200
      X(I)=LOG(XLOWER)+XSTEP*FLOAT(I-1)
      XX(I)=EXP(X(I))
C
C  SWITCH FOR INITIAL CONDITION
C

      F(I,1)=1.
C      F(I,1)=EXP(X(I))

10    CONTINUE

      DO 100 IT=2,61

      DO 110 I=1,200
      DO 110 J=1,200
      B(I,J)=0.
110   CONTINUE

      DO 120 I=1,200
      A(I) =   - F(I,IT-1)/TSTEP
120   CONTINUE

      B(1,1)= -(ALP-BET*X(1))/XSTEP-EXP(X(1))-1./TSTEP
      B(1,2)= (ALP-BET*X(1))/XSTEP
      B(200,199)= -(ALP-BET*X(100))/XSTEP
      B(200,200)= (ALP-BET*X(100))/XSTEP-EXP(X(100))-1./TSTEP

      DO 130 I=2,199
      B(I,I-1)=  SIG**2./(2.*XSTEP**2.)-(ALP-BET*X(I))/(2.*XSTEP)
      B(I,I)=   -SIG**2./XSTEP**2.-EXP(X(I))-1./TSTEP
      B(I,I+1)=  SIG**2./(2.*XSTEP**2.)+(ALP-BET*X(I))/(2.*XSTEP)
130   CONTINUE

      DO 140 I=1,199
      C1(I)=B(I,I+1)
140   CONTINUE

      DO 150 I=1,200
      C2(I)=B(I,I)
150   CONTINUE

      DO 160 I=2,200
      C3(I)=B(I,I-1)
160   CONTINUE

      N=200
      CALL DLSLTR(N,C3,C2,C1,A)

      DO 170 I=1,200
      F(I,IT)=A(I)
170   CONTINUE

100   CONTINUE


1000  FORMAT(3X,50F12.7)

      END

C  SZ END PDE1P ========================================================




C  SZ START PDE2P ========================================================

      SUBROUTINE PDE2P(PARM,XX,G)
      IMPLICIT REAL*8 (A-H,P-Z)

C
C  THIS SUBROUTINE SOLVES THE
C  PDE FOR THE SECOND TYPE OF
C  INITIAL CONDITION AND RETURNS
C  THE MATRIX OF VALUES AS F
C
C  UNDER P MEASURE
C

      DIMENSION X(200),XX(200)
      DIMENSION G(200,61)
      DIMENSION A(200),B(200,200)
      DIMENSION PARM(10)
      DIMENSION C1(200),C2(200),C3(200)

C
C  NEED TO SET LOWER LIMIT
C   AND STEP SIZES
C

      XLOWER = .0001
      XUPPER = .2500

      TSTEP=  1./12.
      XSTEP=  (LOG(XUPPER)-LOG(XLOWER))/200.

      ALP  = PARM(1)
      BET  = PARM(2)
      SIG  = PARM(3)
      ALPH = PARM(4)
      BETH = PARM(5)

      ALP=ALPH
      BET=BETH

      DO 10 I=1,200
      X(I)=LOG(XLOWER)+XSTEP*FLOAT(I-1)
      XX(I)=EXP(X(I))

C
C  SWITCH FOR INITIAL CONDITION
C

C      G(I,1)=1.
      G(I,1)=EXP(X(I))

10    CONTINUE

      DO 100 IT=2,61

      DO 110 I=1,200
      DO 110 J=1,200
      B(I,J)=0.
110   CONTINUE

      DO 120 I=1,200
      A(I) =   - G(I,IT-1)/TSTEP
120   CONTINUE

      B(1,1)= -(ALP-BET*X(1))/XSTEP-EXP(X(1))-1./TSTEP
      B(1,2)= (ALP-BET*X(1))/XSTEP
      B(200,199)= -(ALP-BET*X(100))/XSTEP
      B(200,200)= (ALP-BET*X(100))/XSTEP-EXP(X(100))-1./TSTEP

      DO 130 I=2,199
      B(I,I-1)=  SIG**2./(2.*XSTEP**2.)-(ALP-BET*X(I))/(2.*XSTEP)
      B(I,I)=   -SIG**2./XSTEP**2.-EXP(X(I))-1./TSTEP
      B(I,I+1)=  SIG**2./(2.*XSTEP**2.)+(ALP-BET*X(I))/(2.*XSTEP)
130   CONTINUE

      DO 140 I=1,199
      C1(I)=B(I,I+1)
140   CONTINUE

      DO 150 I=1,200
      C2(I)=B(I,I)
150   CONTINUE

      DO 160 I=2,200
      C3(I)=B(I,I-1)
160   CONTINUE

      N=200
      CALL DLSLTR(N,C3,C2,C1,A)

      DO 170 I=1,200
      G(I,IT)=A(I)
170   CONTINUE

100   CONTINUE


1000  FORMAT(3X,50F12.7)

      END

C  SZ END PDE2P ========================================================







C  SZ START XINTPP ========================================================


      SUBROUTINE XINTPP(CDS3,X,XLAM,
     *                 CDS1T,CDS3T,CDS5T,
     *                 CDS1X,CDS3X,CDS5X)

      IMPLICIT REAL*8 (A-H,P-Z)

C
C    THIS SUBROUTINE INTERPOLATES
C    MATRIX OF CDS PRICES TO MATCH
C    MARKET DATA
C
C    UNDER P MEASURE
C

C  SZ CHANGE DATA, PFILE

      DIMENSION PARM(10)
      DIMENSION DATA(85,15,3)
      DIMENSION PFILE(85,260)
      DIMENSION ICAR(10)
      DIMENSION X(200)
      DIMENSION F(200,61)
      DIMENSION G(200,61)
      DIMENSION CDS1T(200),CDS3T(200),CDS5T(200)

      COMMON ICAR,PARM,DATA,PFILE

      DO 10 I=1,199
      IF((XLAM.GE.X(I)).AND.(XLAM.LE.X(I+1)))THEN

      WEIGHT = (XLAM-X(I))/(X(I+1)-X(I))
C      XLAM   = (1.-WEIGHT)*X(I)+WEIGHT*X(I+1)
      CDS1X  = (1.-WEIGHT)*CDS1T(I)+WEIGHT*CDS1T(I+1)
      CDS3X  = (1.-WEIGHT)*CDS3T(I)+WEIGHT*CDS3T(I+1)
      CDS5X  = (1.-WEIGHT)*CDS5T(I)+WEIGHT*CDS5T(I+1)



C      WRITE(*,1000)CDS3T(I+1),CDS3T(I),LOG(X(I+1)),LOG(X(I))

      GO TO 20
      ENDIF

10    CONTINUE
      XLAM=X(200)
      CDS1X=CDS1T(200)
      CDS3X=CDS3T(200)
      CDS5X=CDS5T(200)


1000  FORMAT(3X,20F16.7)

20    CONTINUE
      RETURN
      END

C  SZ END XINTPP ========================================================

