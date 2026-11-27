C*********************************************************************
C  CIRCULAR-DOMAIN DIRICHLET EIGENMODE RECONSTRUCTION
C
C  COMPLEX FORMULATION:
C      A = [U;L]
C  THE COMPLEX SVD IS PERFORMED THROUGH THE EQUIVALENT REAL MATRIX
C
C      AR = [ REAL(A)  -IMAG(A) ]
C           [ IMAG(A)   REAL(A) ]
C
C  INPUT : K, F15.DAT, 1F80.DAT
C  OUTPUT: xspmoU2.dat,   xspedU2.DAT
C          xspmoU2_2.dat, xspedU2_2.DAT
C          THE SECOND-MODE FILES CONTAIN ZEROS FOR A SINGLE ROOT.
C*********************************************************************
      IMPLICIT REAL*8 (A-H,O-Y)

      PARAMETER(NELMAX=500)
      PARAMETER(MMAX=4*NELMAX,NMAX=2*NELMAX)

      COMPLEX*16 ZVALUEU,ZVALUEL,ZMODE1,ZMODE2
      COMPLEX*16 ZVEC1(NELMAX),ZVEC2(NELMAX)

      COMMON/AA/ NELM,NNODE
      COMMON/NORMAL/ XNORMAL,YNORMAL,THETA
      COMMON/BAR/ XBAR,YBAR,XROT,YROT

      DIMENSION INC(2,NELMAX),ELENGTH(NELMAX)
      DIMENSION THETA(NELMAX),X(NELMAX),Y(NELMAX)
      DIMENSION NINDEX(NELMAX),COLX(NELMAX),COLY(NELMAX)
      DIMENSION XNORMAL(NELMAX),YNORMAL(NELMAX)
      DIMENSION AM(MMAX,NMAX),W1(NMAX),V1(NMAX,NMAX)

      INTEGER ID(4),ID1,ID2,NMODE,IOS,ITEMP
      REAL*8 TOL,RATIO,MODE1,MODE2

      OPEN(15,FILE='F15.DAT',STATUS='OLD',IOSTAT=IOS)
      IF(IOS.NE.0) STOP 'ERROR: CANNOT OPEN F15.DAT.'

      OPEN(80,FILE='1F80.DAT',STATUS='OLD',IOSTAT=IOS)
      IF(IOS.NE.0) STOP 'ERROR: CANNOT OPEN 1F80.DAT.'

      PRINT *,'INPUT EIGENVALUE K'
      READ(5,*) EIG

      CALL COUNT_F15(NNODE,NELM)
      REWIND(15)
      CALL COUNT_INTERIOR(NINTER)
      REWIND(80)

      IF(NNODE.LE.0.OR.NELM.LE.0) STOP 'ERROR: INVALID F15.DAT.'
      IF(NINTER.LE.0) STOP 'ERROR: INVALID 1F80.DAT.'
      IF(NNODE.GT.NELMAX.OR.NELM.GT.NELMAX) THEN
         PRINT *,'ERROR: NELM OR NNODE EXCEEDS',NELMAX
         STOP
      END IF

      PRINT *,'NELM   =',NELM
      PRINT *,'NNODE  =',NNODE
      PRINT *,'NINTER =',NINTER
      PRINT *,'K      =',EIG

      CALL GEOM(X,Y,INC,NINDEX,NELM,NNODE)

C----- COLLOCATION POINTS, ELEMENT LENGTHS AND OUTWARD NORMALS
      DO 10 I=1,NELM
         J=NINDEX(I)
         N1=INC(1,J)
         N2=INC(2,J)

         COLX(J)=0.5D0*(X(N1)+X(N2))
         COLY(J)=0.5D0*(Y(N1)+Y(N2))

         DELTX=X(N2)-X(N1)
         DELTY=Y(N2)-Y(N1)
         ELENGTH(J)=DSQRT(DELTX**2+DELTY**2)

         XNORMAL(J)=DELTY/ELENGTH(J)
         YNORMAL(J)=-DELTX/ELENGTH(J)
         THETA(J)=DATAN2(DELTY,DELTX)*180.0D0/3.141592654D0
   10 CONTINUE

C----- BUILD THE REAL EQUIVALENT OF THE COMPLEX MATRIX [U;L]
      DO 20 ICOL=1,NELM
         DO 20 NEL=1,NELM
            XBAR=COLX(ICOL)-COLX(NEL)
            YBAR=COLY(ICOL)-COLY(NEL)
            THET=THETA(NEL)*3.141592654D0/180.0D0

            XROT=XBAR*DCOS(THET)+YBAR*DSIN(THET)
            YROT=-XBAR*DSIN(THET)+YBAR*DCOS(THET)

            CALL UKERNEL(EIG,NEL,XROT,YROT,ELENGTH,ZVALUEU)
            CALL LKERNEL(EIG,ICOL,NEL,XROT,YROT,ELENGTH,
     *                   ZVALUEL)

C----------- U BLOCK
            AM(ICOL,NEL)=DREAL(ZVALUEU)
            AM(ICOL,NEL+NELM)=-DIMAG(ZVALUEU)
            AM(ICOL+2*NELM,NEL)=DIMAG(ZVALUEU)
            AM(ICOL+2*NELM,NEL+NELM)=DREAL(ZVALUEU)

C----------- L BLOCK
            IROW=ICOL+NELM
            AM(IROW,NEL)=DREAL(ZVALUEL)
            AM(IROW,NEL+NELM)=-DIMAG(ZVALUEL)
            AM(IROW+2*NELM,NEL)=DIMAG(ZVALUEL)
            AM(IROW+2*NELM,NEL+NELM)=DREAL(ZVALUEL)
   20 CONTINUE

      CALL SVDCMP(AM,4*NELM,2*NELM,MMAX,NMAX,W1,V1)
      CALL FOUR_SMALLEST(W1,2*NELM,ID)

C----- EACH COMPLEX SINGULAR VALUE APPEARS TWICE IN THE REAL MATRIX
      ID1=ID(1)
      ID2=ID(3)

      TOL=1.0D-3
      RATIO=DABS(W1(ID2)-W1(ID1))
     *      /DMAX1(DABS(W1(ID1)),1.0D-30)
      NMODE=1
      IF(RATIO.LT.TOL) NMODE=2

      CALL MAKE_COMPLEX_VECTOR(V1,NMAX,NELM,ID1,ZVEC1)
      IF(NMODE.EQ.2) THEN
         CALL MAKE_COMPLEX_VECTOR(V1,NMAX,NELM,ID2,ZVEC2)
      ELSE
         DO 25 I=1,NELM
            ZVEC2(I)=(0.0D0,0.0D0)
   25    CONTINUE
      END IF

      PRINT *,'COMPLEX SV1 =',W1(ID1)
      PRINT *,'COMPLEX SV2 =',W1(ID2)
      PRINT *,'NMODE       =',NMODE

      OPEN(81,FILE='xspmoU2.dat',STATUS='REPLACE')
      OPEN(300,FILE='xspedU2.DAT',STATUS='REPLACE')
      OPEN(812,FILE='xspmoU2_2.dat',STATUS='REPLACE')
      OPEN(302,FILE='xspedU2_2.DAT',STATUS='REPLACE')

C----- REAL REPRESENTATIVES OF THE COMPLEX BOUNDARY DENSITIES
      DO 30 I=1,NELM
         WRITE(300,102) I,DREAL(ZVEC1(I))
         IF(NMODE.EQ.2) THEN
            WRITE(302,102) I,DREAL(ZVEC2(I))
         ELSE
            WRITE(302,102) I,0.0D0
         END IF
   30 CONTINUE

C----- RECONSTRUCT THE COMPLEX INTERIOR MODES, THEN TAKE REAL PART
      DO 40 I=1,NINTER
         READ(80,100) ITEMP,XP,YP
         ZMODE1=(0.0D0,0.0D0)
         ZMODE2=(0.0D0,0.0D0)

         DO 41 NEL=1,NELM
            XBAR=XP-COLX(NEL)
            YBAR=YP-COLY(NEL)
            THET=THETA(NEL)*3.141592654D0/180.0D0

            XROT=XBAR*DCOS(THET)+YBAR*DSIN(THET)
            YROT=-XBAR*DSIN(THET)+YBAR*DCOS(THET)

            CALL UKERNEL(EIG,NEL,XROT,YROT,ELENGTH,ZVALUEU)
            ZMODE1=ZMODE1-ZVALUEU*ZVEC1(NEL)
            IF(NMODE.EQ.2) THEN
               ZMODE2=ZMODE2-ZVALUEU*ZVEC2(NEL)
            END IF
   41    CONTINUE

         ZMODE1=ZMODE1/(2.0D0*3.141592654D0)
         MODE1=DREAL(ZMODE1)
         WRITE(81,101) XP,YP,MODE1

         IF(NMODE.EQ.2) THEN
            ZMODE2=ZMODE2/(2.0D0*3.141592654D0)
            MODE2=DREAL(ZMODE2)
         ELSE
            MODE2=0.0D0
         END IF
         WRITE(812,101) XP,YP,MODE2
   40 CONTINUE

      PRINT *,'COMPLEX MODE RECONSTRUCTION FINISHED.'

      CLOSE(15)
      CLOSE(80)
      CLOSE(81)
      CLOSE(300)
      CLOSE(812)
      CLOSE(302)

  100 FORMAT(I6,2E13.5)
  101 FORMAT(3F13.5)
  102 FORMAT(I6,F13.5)
      STOP
      END

C*********************************************************************
C  COUNT NODES AND ELEMENTS IN F15.DAT
C*********************************************************************
      SUBROUTINE COUNT_F15(NNODE,NELM)
      IMPLICIT REAL*8 (A-H,O-Y)
      CHARACTER*80 LINE
      INTEGER IOS,IVAL,NLINE

      NNODE=0
      NELM=0
      REWIND(15)

   10 READ(15,'(A80)',IOSTAT=IOS) LINE
      IF(IOS.NE.0) RETURN
      IF(LINE(4:6).NE.' 15') GO TO 10

   20 READ(15,'(A80)',IOSTAT=IOS) LINE
      IF(IOS.NE.0) RETURN
      READ(LINE,*,IOSTAT=IOS) IVAL
      IF(IOS.NE.0) GO TO 20
      IF(IVAL.EQ.-1) GO TO 30
      NNODE=NNODE+1
      GO TO 20

   30 READ(15,'(A80)',IOSTAT=IOS) LINE
      IF(IOS.NE.0) RETURN
      IF(LINE(4:6).NE.' 71') GO TO 30

      NLINE=0
   40 READ(15,'(A80)',IOSTAT=IOS) LINE
      IF(IOS.NE.0) GO TO 50
      READ(LINE,*,IOSTAT=IOS) IVAL
      IF(IOS.NE.0) GO TO 40
      IF(IVAL.EQ.-1) GO TO 50
      NLINE=NLINE+1
      GO TO 40

   50 NELM=NLINE/2
      RETURN
      END

C*********************************************************************
C  COUNT INTERIOR POINTS IN 1F80.DAT
C*********************************************************************
      SUBROUTINE COUNT_INTERIOR(NINTER)
      IMPLICIT REAL*8 (A-H,O-Y)
      INTEGER IOS,ITEMP

      NINTER=0
   10 READ(80,*,IOSTAT=IOS) ITEMP,XP,YP
      IF(IOS.NE.0) RETURN
      NINTER=NINTER+1
      GO TO 10
      END

C*********************************************************************
C  INDICES OF THE FOUR SMALLEST SINGULAR VALUES
C*********************************************************************
      SUBROUTINE FOUR_SMALLEST(W,N,ID)
      IMPLICIT REAL*8 (A-H,O-Y)
      DIMENSION W(N)
      INTEGER ID(4),ITMP

      DO 10 I=1,4
         ID(I)=I
   10 CONTINUE

      DO 30 I=1,4
         DO 20 J=I+1,4
            IF(DABS(W(ID(J))).LT.DABS(W(ID(I)))) THEN
               ITMP=ID(I)
               ID(I)=ID(J)
               ID(J)=ITMP
            END IF
   20    CONTINUE
   30 CONTINUE

      DO 50 I=5,N
         IF(DABS(W(I)).LT.DABS(W(ID(4)))) THEN
            ID(4)=I
            DO 40 J=4,2,-1
               IF(DABS(W(ID(J))).LT.DABS(W(ID(J-1)))) THEN
                  ITMP=ID(J)
                  ID(J)=ID(J-1)
                  ID(J-1)=ITMP
               END IF
   40       CONTINUE
         END IF
   50 CONTINUE
      RETURN
      END

C*********************************************************************
C  CONVERT ONE REALIFIED RIGHT SINGULAR VECTOR TO A COMPLEX VECTOR
C  V(:,ID) = [REAL(ZVEC); IMAG(ZVEC)]
C  THE GLOBAL PHASE IS ADJUSTED SO THE LARGEST COMPONENT IS REAL.
C*********************************************************************
      SUBROUTINE MAKE_COMPLEX_VECTOR(V,LD,N,ID,ZVEC)
      IMPLICIT REAL*8 (A-H,O-Y)
      INTEGER LD,N,ID,IMAX
      REAL*8 V(LD,LD),VMAX
      COMPLEX*16 ZVEC(N),ZPHASE

      VMAX=0.0D0
      IMAX=1
      DO 10 I=1,N
         ZVEC(I)=DCMPLX(V(I,ID),V(I+N,ID))
         IF(CDABS(ZVEC(I)).GT.VMAX) THEN
            VMAX=CDABS(ZVEC(I))
            IMAX=I
         END IF
   10 CONTINUE

      IF(VMAX.GT.1.0D-30) THEN
         ZPHASE=DCONJG(ZVEC(IMAX))/VMAX
         DO 20 I=1,N
            ZVEC(I)=ZVEC(I)*ZPHASE
   20    CONTINUE
      END IF
      RETURN
      END

      SUBROUTINE UKERNEL(EIG,NEL,XROT,YROT,ELENGTH,ZVALUEU)
       IMPLICIT REAL*8 (A-H,O-Y)
           COMPLEX*16 ZU,ZU0,ZVALUEU,ZV
           DIMENSION XW(20),W(20),ELENGTH(500)
           save
C*****************************************************
C         FOR SINGULAR ELEMENT                       *
C*****************************************************
        IF(DABS(XROT).LE.0.00000001.AND.DABS(YROT).LE.0.00000001) THEN
              R=EIG*0.5d0*ELENGTH(NEL)
             CALL BESSH0(R,HJ0,HY0)
            ZU0=HJ0+(0.d0,1.d0)*HY0
           ZV=ZU0*ELENGTH(NEL)
          CALL GRULE(20,XW,W)
          CALL ARANGE(20,XW,W)
           ZVALUEU=(0.d0,0.d0)
          DO 555 I=1,20
              X=DABS(0.5D0*ELENGTH(NEL)*XW(I))
              R=EIG*X
             CALL BESSH1(R,HJ1,HY1)
            ZU=HJ1+(0.d0,1.d0)*HY1
            ZVALUEU=ZVALUEU+ZU*X*W(I)
 555      CONTINUE
            ZVALUEU=ZVALUEU*EIG*0.5d0*ELENGTH(NEL)

         ZVALUEU=(ZV+ZVALUEU)*(-3.141592654d0/2.d0)*(0.d0,1.d0)
C*****************************************************
C         FOR REGULAR ELEMENT                       *
C*****************************************************
        ELSE
          CALL GRULE(20,XW,W)
          CALL ARANGE(20,XW,W)
           ZVALUEU=(0.,0.)
          DO 655 I=1,20
              R=EIG*DSQRT(YROT**2+(0.5D0*ELENGTH(NEL)*XW(I)-XROT)**2)
             CALL BESSH0(R,HJ0,HY0)
            ZU=HJ0+(0.D0,1.D0)*HY0
            ZVALUEU=ZVALUEU+ZU*W(I)*0.5D0*ELENGTH(NEL)
 655      CONTINUE
         ZVALUEU=ZVALUEU*(-3.141592654D0/2.D0)*(0.D0,1.D0)
        ENDIF
        RETURN
        END

CC**************************************************************
C    TKERNEL SUBROUTINE CAN CALCULATE THE T(s,x) INTEGRATION  *
C      INCLUDING SINGULAR AND NONSINGULAR VALUE               *
C                                                             *
C**************************************************************
      SUBROUTINE LKERNEL(EIG,ICOL,NEL,XROT,YROT,ELENGTH,ZVALUEL)
       IMPLICIT REAL*8 (A-H,O-Y)
       COMPLEX*16 ZU,ZVALUEL
      COMMON /NORMAL/ XNORMAL,YNORMAL,THETA
      COMMON/AA/ NELM,NNODE
      DIMENSION XW(20),W(20),ELENGTH(500)
     *          ,XNORMAL(500),YNORMAL(500)
     *          ,THETA(500)
         XN=XNORMAL(ICOL)
         XS=XNORMAL(NEL)
         YN=YNORMAL(ICOL)
         YS=YNORMAL(NEL)
         CDAR=THETA(NEL)*3.141592654D0/180.D0
         DOT=XN*XS+YN*YS
         CROSS=-1.D0*XN*YS+YN*XS
C*****************************************************
C         FOR SINGULAR ELEMENT                       *
C*****************************************************
        IF(ABS(XROT).LE.0.00000001.AND.ABS(YROT).LE.0.00000001) THEN
           ZVALUEL=(3.141592654D0,0.D0)
        IF(ABS(XROT).LE.0.00000001.AND.ICOL.NE.NEL ) ZVALUEL=(-3.14159
     C 2654D0,0.D0)
       ELSE
C*****************************************************
C         FOR REGULAR ELEMENT                       *
C*****************************************************
          CALL GRULE(20,XW,W)
          CALL ARANGE(20,XW,W)
           ZVALUEL=(0.,0.)
          DO 755 I=1,20
             X=DSQRT(YROT**2+(0.5D0*ELENGTH(NEL)*XW(I)-XROT)**2)
             R=EIG*X
             CALL BESSH1(R,HJ1,HY1)
          ZU=HJ1+(0.D0,1.D0)*HY1
          ZVALUEL=ZVALUEL+ZU*((XROT-0.5D0*ELENGTH(NEL)*XW(I))*CROSS-
     *          YROT*DOT)/X*W(I)
 755      CONTINUE
          ZVALUEL=ZVALUEL*0.5D0*ELENGTH(NEL)
         ZVALUEL=ZVALUEL*(3.141592654D0/2.D0)*(0.D0,1.D0)*EIG
       ENDIF
       RETURN
       END
C****************************************************

CC*************************************************************
C    MKERNEL SUBROUTINE CAN CALCULATE THE M(s,x) INTEGRATION  *
C      INCLUDING SINGULAR AND NONSINGULAR VALUE               *
C         BY GAUSSIAN INTEGRATION FORM SOLUTION               *
C**************************************************************
C         SUBROUTINE GRULE                             *
C   THIS IS A SUBROUTINE TO COMPUTE THE (N+1)/2        *
C   NONNEGATIVE AB-SCISSAS X(I),W(I) OF THE N-PT       *
C   GAUSSIAN-LENGENDRE INTEGRATION RULE NORMALIZED     *
C   TO THE INTERVAL [-1,1]                             *
C*******************************************************
          SUBROUTINE GRULE(N,X,W)
      IMPLICIT REAL*8 (A-H,O-Y)
          DIMENSION X(N),W(N)
          M=(N+1)/2
          E1=N*(N+1)
          DO 1 I=1,M
            T=(4.0*I-1.0)*3.1415926536/(4*N+2)
            X0=(1.-(1.-1./N)/(8.*N*N))*COS(T)
            PKM1=1.0
            PK=X0
            DO 3 K=2,N
              T1=PK*X0
              PKP1=T1-PKM1-(T1-PKM1)/K+T1
              PKM1=PK
              PK=PKP1
3           CONTINUE
            DEN=1.-X0*X0
            D1=N*(PKM1-X0*PK)
            DPN=D1/DEN
            D2PN=(2.*X0*DPN-E1*PK)/DEN
            D3PN=(4.*X0*D2PN+(2.-E1)*DPN)/DEN
            D4PN=(6.*X0*D3PN+(6.-E1)*D2PN)/DEN
            U=PK/DPN
            V=D2PN/DPN
            H=-U*(1.+.5*U*(V+U*(V*V-U*D3PN/(3.*DPN))))
            P=PK+H*(DPN+0.5*H*(D2PN+H/3.*(D3PN+.25*H*D4PN)))
            DP=DPN+H*(D2PN+.5*H*(D3PN+H*D4PN/3.))
            H=H-P/DP
            X(I)=X0+H
           FX=D1-H*E1*(PK+.5*H*(DPN+H/3.*(D2PN+.25*H*(D3PN+.2*H*D4PN))))
           W(I)=2.*(1.-X(I)*X(I))/(FX*FX)
1        CONTINUE
         IF(M+M.GT.N)X(M)=0.
         RETURN
         END
C
C
C
C***********************************************
C          SUBROUTINE ARANGE                   *
C  THE SUBROUTINE  IS  TO ARANGE THE GAUSSIAN  *
C  QUA COORD. AND WEIGHTING FUNCTION           *
C***********************************************
      SUBROUTINE ARANGE(N,X,W)
       IMPLICIT REAL*8 (A-H,O-Y)
      REAL*8 X(N),W(N)
      M=N/2
      DO 1 I=1,M
        X(N-I+1)=-X(I)
        W(N-I+1)=W(I)
1     CONTINUE
      RETURN
      END

C**===========================**C
C**   N' ORDER HANK FUNCTION  **C
C**===========================**C

      SUBROUTINE BESSH0(X,HJ0,HY0)
       IMPLICIT REAL*8 (A-H,O-Y)
C**   0' ORDER HANK FUNCTION -C.W.LIN
      IF(X-3.0) 1,1,3
    1 XX=X*X/9.00000000
      HJ0=  1.0000000-(2.2499997-(1.2656208-(0.3163866-(0.0444479
     *               -(0.0039444-0.0002100*XX)*XX)*XX)*XX)*XX)*XX
      HY0=  .63661977*DLOG(0.5*X)*HJ0
     *     +.36746691+(.60559366-(.74350384-(.25300117-(.04261214
     *               -(.00427916-.00024846*XX)*XX)*XX)*XX)*XX)*XX
      RETURN
    3 XX=3.00000000/X
      F0=   .79788456-(.00000077+(.00552740+(.00009512-(.00137237
     *               -(.00072805-.00014476*XX)*XX)*XX)*XX)*XX)*XX
      G0= X-.78539816-(.04166397+(.00003954-(.00262573-(.00054125
     *               +(.00029333-.00013558*XX)*XX)*XX)*XX)*XX)*XX
      XF=F0/DSQRT(X)
      HJ0=XF*DCOS(G0)
      HY0=XF*DSIN(G0)
      RETURN
      END
C**=================================**C
C**     HANKEL FUNCTION             **C
C**        ORDER:1                  **C
C**=================================**C
      SUBROUTINE BESSH1(X,HJ1,HY1)
       IMPLICIT REAL*8 (A-H,O-Y)
C**   1' ORDER HANK FUNCTION -C.W.LIN
      IF(X-3.0) 1,1,3
    1 XX=X*X/9.00000000
      HJ1= (.50000000-(.56249985-(.21093573-(.03954289-(.00443319
     *               -(.00031761-.00001109*XX)*XX)*XX)*XX)*XX)*XX)*X
      HY1=  .63661977*DLOG(0.5*X)*HJ1
     *    -(0.6366198-(0.2212091+(2.1682709-(1.3164827-(0.3123951
     *               -(0.0400967-0.0027873*XX)*XX)*XX)*XX)*XX)*XX)/X
      RETURN
    3 XX=3.00000000/X
      F1=   .79788456+(.00000156+(.01659667+(.00017105-(.00249511
     *               -(.00113653-.00020033*XX)*XX)*XX)*XX)*XX)*XX
      G1=X-2.35619449+(.12499612+(.00005650-(.00637879-(.00074348
     *               +(.00079824-.00029166*XX)*XX)*XX)*XX)*XX)*XX
      XF=F1/DSQRT(X)
      HJ1=XF*DCOS(G1)
      HY1=XF*DSIN(G1)
      RETURN
      END
      SUBROUTINE GEOM(X,Y,INC,NINDEX,NELM,NNODE)
       IMPLICIT REAL*8 (A-H,O-Y)
      DIMENSION X(500),Y(500),NINDEX(500)
      DIMENSION INC(2,500)
      CHARACTER*80 TITLE0,TITLE1
121   FORMAT(36A)
101   READ(15,501)TITLE0
      IF(TITLE0(4:6).NE.' 15')GO TO 101
C*****READ THE CORNER POINT COORDINATES X,Y,Z **********************
      DO 7001 K=1,NNODE
        READ(15,100,END=7001) M1,M2,M3,M4,PP,Q,R
        X(M1)=PP
        Y(M1)=Q
100     FORMAT(4I10,3E13.5)
501     FORMAT(80A)
7001  CONTINUE
202   READ(15,501)TITLE1
      IF(TITLE1(4:6).NE.' 71')GO TO 202
C*****READ THE INCIDENCE OF EACH ELEMENT ***************************
      DO 8001 J=1,NELM
        READ(15,800) II,JJ,KK,LL,MM,NN,IL
        NINDEX(J)=II
800     FORMAT(7I10)
        READ(15,200,END=8001) NINC1,NINC2
        INC(1,II)=NINC1
        INC(2,II)=NINC2
200     FORMAT(2I10)
8001  CONTINUE
      RETURN
      END


      SUBROUTINE SVDCMP(AA,M,N,MP,NP,W,V)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NMAX=1000)
      DIMENSION AA(MP,NP),W(NP),V(NP,NP),RV1(NMAX)
      G=0.0D0
      SCALE=0.0D0
      ANORM=0.0D0
      DO 25 I=1,N
        L=I+1
        RV1(I)=SCALE*G
        G=0.0D0
        S=0.0D0
        SCALE=0.0D0
        IF (I.LE.M) THEN
          DO 11 K=I,M
            SCALE=SCALE+DABS(AA(K,I))
 11        CONTINUE
          IF (SCALE.NE.0.0D0) THEN
            DO 12 K=I,M
              AA(K,I)=AA(K,I)/SCALE
              S=S+AA(K,I)*AA(K,I)
 12          CONTINUE
            F=AA(I,I)
            G=-DSIGN(DSQRT(S),F)
            H=F*G-S
            AA(I,I)=F-G
            IF (I.NE.N) THEN
              DO 15 J=L,N
                S=0.0D0
                DO 13 K=I,M
                  S=S+AA(K,I)*AA(K,J)
 13              CONTINUE
                F=S/H
                DO 14 K=I,M
                  AA(K,J)=AA(K,J)+F*AA(K,I)
 14              CONTINUE
 15            CONTINUE
            ENDIF
            DO 16 K= I,M
              AA(K,I)=SCALE*AA(K,I)
 16          CONTINUE
          ENDIF
        ENDIF
        W(I)=SCALE *G
        G=0.0D0
        S=0.0D0
        SCALE=0.0D0
        IF ((I.LE.M).AND.(I.NE.N)) THEN
          DO 17 K=L,N
            SCALE=SCALE+DABS(AA(I,K))
 17        CONTINUE
          IF (SCALE.NE.0.0D0) THEN
            DO 18 K=L,N
              AA(I,K)=AA(I,K)/SCALE
              S=S+AA(I,K)*AA(I,K)
 18          CONTINUE
            F=AA(I,L)
            G=-DSIGN(DSQRT(S),F)
            H=F*G-S
            AA(I,L)=F-G
            DO 19 K=L,N
              RV1(K)=AA(I,K)/H
 19          CONTINUE
            IF (I.NE.M) THEN
              DO 23 J=L,M
                S=0.0D0
                DO 21 K=L,N
                  S=S+AA(J,K)*AA(I,K)
 21              CONTINUE
                DO 22 K=L,N
                  AA(J,K)=AA(J,K)+S*RV1(K)
 22              CONTINUE
 23            CONTINUE
            ENDIF
            DO 24 K=L,N
              AA(I,K)=SCALE*AA(I,K)
 24          CONTINUE
          ENDIF
        ENDIF
        ANORM=DMAX1(ANORM,(DABS(W(I))+DABS(RV1(I))))
 25    CONTINUE
      DO 32 I=N,1,-1
        IF (I.LT.N) THEN
          IF (G.NE.0.0D0) THEN
            DO 26 J=L,N
              V(J,I)=(AA(I,J)/AA(I,L))/G
 26          CONTINUE
            DO 29 J=L,N
              S=0.0D0
              DO 27 K=L,N
                S=S+AA(I,K)*V(K,J)
 27            CONTINUE
              DO 28 K=L,N
                V(K,J)=V(K,J)+S*V(K,I)
 28            CONTINUE
 29          CONTINUE
          ENDIF
          DO 31 J=L,N
            V(I,J)=0.0D0
            V(J,I)=0.0D0
 31        CONTINUE
        ENDIF
        V(I,I)=1.0D0
        G=RV1(I)
        L=I
 32    CONTINUE
      DO 39 I=N,1,-1
        L=I+1
        G=W(I)
        IF (I.LT.N) THEN
          DO 33 J=L,N
            AA(I,J)=0.0D0
 33        CONTINUE
        ENDIF
        IF (G.NE.0.0D0) THEN
          G=1.0D0/G
          IF (I.NE.N) THEN
            DO 36 J=L,N
              S=0.0D0
              DO 34 K=L,M
                S=S+AA(K,I)*AA(K,J)
 34            CONTINUE
              F=(S/AA(I,I))*G
              DO 35 K=I,M
                AA(K,J)=AA(K,J)+F*AA(K,I)
 35            CONTINUE
 36          CONTINUE
          ENDIF
          DO 37 J=I,M
            AA(J,I)=AA(J,I)*G
 37        CONTINUE
        ELSE
          DO 38 J= I,M
            AA(J,I)=0.0D0
 38        CONTINUE
        ENDIF
        AA(I,I)=AA(I,I)+1.0D0
 39    CONTINUE
      DO 49 K=N,1,-1
        DO 48 ITS=1,30
          DO 41 L=K,1,-1
            NM=L-1
            IF ((DABS(RV1(L))+ANORM).EQ.ANORM)  GO TO 2
            IF ((DABS(W(NM))+ANORM).EQ.ANORM)  GO TO 1
 41        CONTINUE
 1         C=0.0D0
          S=1.0D0
          DO 43 I=L,K
            F=S*RV1(I)
            IF ((DABS(F)+ANORM).NE.ANORM) THEN
              G=W(I)
              H=DSQRT(F*F+G*G)
              W(I)=H
              H=1.0D0/H
              C= (G*H)
              S=-(F*H)
              DO 42 J=1,M
                Y=AA(J,NM)
                Z=AA(J,I)
                AA(J,NM)=(Y*C)+(Z*S)
                AA(J,I)=-(Y*S)+(Z*C)
 42            CONTINUE
            ENDIF
 43        CONTINUE
 2         Z=W(K)
          IF (L.EQ.K) THEN
            IF (Z.LT.0.0D0) THEN
              W(K)=-Z
              DO 44 J=1,N
                V(J,K)=-V(J,K)
 44            CONTINUE
            ENDIF
            GO TO 3
          ENDIF
          IF (ITS.EQ.30) PAUSE 'No convergence in 30 iterations'
          X=W(L)
          NM=K-1
          Y=W(NM)
          G=RV1(NM)
          H=RV1(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0D0*H*Y)
          G=DSQRT(F*F+1.0D0)
          F=((X-Z)*(X+Z)+H*((Y/(F+DSIGN(G,F)))-H))/X
          C=1.0D0
          S=1.0D0
          DO 47 J=L,NM
            I=J+1
            G=RV1(I)
            Y=W(I)
            H=S*G
            G=C*G
            Z=DSQRT(F*F+H*H)
            RV1(J)=Z
            C=F/Z
            S=H/Z
            F= (X*C)+(G*S)
            G=-(X*S)+(G*C)
            H=Y*S
            Y=Y*C
            DO 45 NM=1,N
              X=V(NM,J)
              Z=V(NM,I)
              V(NM,J)= (X*C)+(Z*S)
              V(NM,I)=-(X*S)+(Z*C)
 45          CONTINUE
            Z=DSQRT(F*F+H*H)
            W(J)=Z
            IF (Z.NE.0.0D0) THEN
              Z=1.0D0/Z
              C=F*Z
              S=H*Z
            ENDIF
            F= (C*G)+(S*Y)
            X=-(S*G)+(C*Y)
            DO 46 NM=1,M
              Y=AA(NM,J)
              Z=AA(NM,I)
              AA(NM,J)= (Y*C)+(Z*S)
              AA(NM,I)=-(Y*S)+(Z*C)
 46          CONTINUE
 47        CONTINUE
          RV1(L)=0.0D0
          RV1(K)=F
          W(K)=X
 48      CONTINUE
 3       CONTINUE
 49    CONTINUE
*
*       --- end of the subroutine ---
*
      RETURN
      END

C     *********************************************************
