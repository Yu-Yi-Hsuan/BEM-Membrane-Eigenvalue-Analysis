C=====================================================================
C  (N,K)-GON: COMPLEX(U) ONLY
C  OUTPUT: K, MINIMUM SINGULAR VALUE OF COMPLEX(U)
C=====================================================================
      IMPLICIT REAL*8 (A-H,O-Y)
      COMPLEX*16 ZVALUEU

C----- USER SETTINGS
      PARAMETER(NELMAX=500)
      PARAMETER(MMAX=2*NELMAX,NMAX=2*NELMAX)
      PARAMETER(EIGMIN=0.001D0,EIGMAX=6.0D0,DEIG=0.001D0)

      DIMENSION INC(2,NELMAX),ELENGTH(NELMAX)
      DIMENSION THETA(NELMAX),COLX(NELMAX),COLY(NELMAX)
      DIMENSION X(NELMAX),Y(NELMAX),NINDEX(NELMAX)
      COMPLEX*16 ZU(NELMAX,NELMAX)
      DIMENSION AM(MMAX,NMAX),W(NMAX),V(NMAX,NMAX)

      OPEN(15,FILE='F15.DAT',STATUS='OLD')
      OPEN(31,FILE='Complex_U.DAT',STATUS='REPLACE')

      PI=4.0D0*DATAN(1.0D0)

C----- READ F15.DAT
      CALL COUNT_F15(NNODE,NELM)
      IF(NNODE.GT.NELMAX .OR. NELM.GT.NELMAX) THEN
         WRITE(*,*) 'ERROR: F15 SIZE EXCEEDS NELMAX =',NELMAX
         WRITE(*,*) 'NNODE =',NNODE,'  NELM =',NELM
         STOP
      ENDIF
      REWIND(15)

      WRITE(*,*) '(N,K)-GON MODEL: COMPLEX(U) ONLY'
      WRITE(*,*) 'NELM  =',NELM
      WRITE(*,*) 'NNODE =',NNODE

      CALL GEOM(X,Y,INC,NINDEX,NELM,NNODE)

C----- FORM BOUNDARY GEOMETRY
      DO 100 I=1,NELM
         J=NINDEX(I)

         COLX(J)=0.5D0*(X(INC(1,J))+X(INC(2,J)))
         COLY(J)=0.5D0*(Y(INC(1,J))+Y(INC(2,J)))

         DELTX=X(INC(2,J))-X(INC(1,J))
         DELTY=Y(INC(2,J))-Y(INC(1,J))
         ELENGTH(J)=DSQRT(DELTX**2+DELTY**2)

         IF(DABS(DELTY).LE.1.0D-8) GO TO 20
         ARCCOT=90.0D0-DATAN(DELTX/DELTY)*180.0D0/PI
         THETA(J)=ARCCOT
         IF(DELTY.LT.0.0D0) THETA(J)=ARCCOT+180.0D0
         GO TO 30

   20    IF(DELTX.GE.0.0D0) THEN
            THETA(J)=0.0D0
         ELSE
            THETA(J)=180.0D0
         ENDIF

   30    CONTINUE
  100 CONTINUE

      NSTEP=INT((EIGMAX-EIGMIN)/DEIG)
      WRITE(*,*) 'START SCANNING K ...'

C----- SCAN WAVE NUMBER K
      DO 500 ISTEP=0,NSTEP
         EIG=EIGMIN+DBLE(ISTEP)*DEIG

C-------- BUILD COMPLEX U
         DO 200 ICOL=1,NELM
            DO 210 NEL=1,NELM
               XBAR=COLX(ICOL)-COLX(NEL)
               YBAR=COLY(ICOL)-COLY(NEL)

               THET=THETA(NEL)*PI/180.0D0
               XROT=XBAR*DCOS(THET)+YBAR*DSIN(THET)
               YROT=-XBAR*DSIN(THET)+YBAR*DCOS(THET)

               CALL UKERNEL(EIG,NEL,XROT,YROT,ELENGTH,
     *                      ZVALUEU)

               ZU(ICOL,NEL)=ZVALUEU
  210       CONTINUE
  200    CONTINUE

C-------- 1. COMPLEX U
         MROW=2*NELM
         NCOL=2*NELM

         DO 300 I=1,NELM
            DO 310 J=1,NELM
               AM(I,J)=DBLE(ZU(I,J))
               AM(I,J+NELM)=-DIMAG(ZU(I,J))
               AM(I+NELM,J)=DIMAG(ZU(I,J))
               AM(I+NELM,J+NELM)=DBLE(ZU(I,J))
  310       CONTINUE
  300    CONTINUE

         CALL MIN_SV(AM,MROW,NCOL,MMAX,NMAX,W,V,SIGU)
         WRITE(31,*) EIG,SIGU

C-------- DISPLAY PROGRESS EVERY 10 STEPS
         IF(MOD(ISTEP,10).EQ.0) THEN
            WRITE(*,900) EIG,SIGU
         ENDIF

  500 CONTINUE

  900 FORMAT(' K=',F8.3,'  U=',1PE10.3)

      CLOSE(15)
      CLOSE(31)

      WRITE(*,*) 'SCAN FINISHED.'
      WRITE(*,*) 'OUTPUT FILES:'
      WRITE(*,*) '  Complex_U.DAT'
      STOP
      END

C----- RETURN MINIMUM SINGULAR VALUE
      SUBROUTINE MIN_SV(A,M,N,MMAX,NMAX,W,V,SIG)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(MMAX,NMAX),W(NMAX),V(NMAX,NMAX)

      CALL SVDCMP(A,M,N,MMAX,NMAX,W,V)
      SIG=DABS(W(1))

      DO 10 I=2,N
         SIG=DMIN1(SIG,DABS(W(I)))
   10 CONTINUE

      RETURN
      END

C----- COUNT NODES AND ELEMENTS IN F15.DAT
      SUBROUTINE COUNT_F15(NNODE,NELM)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*200 LINE,TMP
      INTEGER IOS,STATE,EXPECT2

      NNODE=0
      NELM=0
      STATE=0
      EXPECT2=0
      REWIND(15)

   10 READ(15,'(A)',IOSTAT=IOS) LINE
      IF(IOS.NE.0) GO TO 90

C----- SECTION MARKERS MUST OCCUPY THE LINE ALONE
      TMP=ADJUSTL(LINE)
      IF(STATE.EQ.0) THEN
         IF(TRIM(TMP).EQ.'15') STATE=1
         GO TO 10
      ENDIF

      IF(STATE.EQ.1) THEN
         READ(LINE,*,IOSTAT=IOS) I1,I2,I3,I4,A1,A2,A3
         IF(IOS.EQ.0) THEN
            NNODE=NNODE+1
         ELSE
            IF(TRIM(TMP).EQ.'71') THEN
               STATE=2
               EXPECT2=0
            ENDIF
         ENDIF
         GO TO 10
      ENDIF

      IF(STATE.EQ.2) THEN
         IF(EXPECT2.EQ.0) THEN
            READ(LINE,*,IOSTAT=IOS) I1,I2,I3,I4,I5,I6,I7
            IF(IOS.EQ.0) THEN
               NELM=NELM+1
               EXPECT2=1
            ENDIF
         ELSE
            READ(LINE,*,IOSTAT=IOS) I1,I2
            IF(IOS.EQ.0) EXPECT2=0
         ENDIF
         GO TO 10
      ENDIF

   90 IF(NNODE.LE.0 .OR. NELM.LE.0) THEN
         WRITE(*,*) 'ERROR: CANNOT COUNT F15.DAT.'
         STOP
      ENDIF
      REWIND(15)
      RETURN
      END

C----- READ F15.DAT
      SUBROUTINE GEOM(X,Y,INC,NINDEX,NELM,NNODE)
      IMPLICIT REAL*8 (A-H,O-Y)
      DIMENSION X(*),Y(*),NINDEX(*)
      DIMENSION INC(2,*)

      READ(15,*) IDUMMY
      READ(15,*) ISECTION

      IF(ISECTION.NE.15) THEN
         WRITE(*,*) 'ERROR: F15 NODE SECTION NOT FOUND.'
         STOP
      ENDIF

      DO 100 K=1,NNODE
         READ(15,10) M1,M2,M3,M4,PP,Q,R
         X(M1)=PP
         Y(M1)=Q
  100 CONTINUE

      READ(15,*) IDUMMY
      READ(15,*) IDUMMY
      READ(15,*) ISECTION

      IF(ISECTION.NE.71) THEN
         WRITE(*,*) 'ERROR: F15 ELEMENT SECTION NOT FOUND.'
         STOP
      ENDIF

      DO 200 J=1,NELM
         READ(15,20) II,JJ,KK,LL,MM,NN,IL
         READ(15,30) NINC1,NINC2

         IF(II.LT.1 .OR. II.GT.NELM) THEN
            WRITE(*,*) 'ERROR: INVALID ELEMENT ID =',II
            STOP
         ENDIF

         NINDEX(J)=II
         INC(1,II)=NINC1
         INC(2,II)=NINC2
  200 CONTINUE

   10 FORMAT(4I10,3E13.5)
   20 FORMAT(7I10)
   30 FORMAT(2I10)

      WRITE(*,*) 'F15.DAT READ COMPLETE.'
      RETURN
      END

      SUBROUTINE UKERNEL(EIG,NEL,XROT,YROT,ELENGTH,ZVALUEU)
       IMPLICIT REAL*8 (A-H,O-Y)
           COMPLEX*16 ZU,ZU0,ZVALUEU,ZV
           DIMENSION XW(20),W(20),ELENGTH(*)
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
           ZVALUEU=(0.D0,0.D0)
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

      SUBROUTINE SVDCMP(AA,M,N,MP,NP,W,V)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AA(MP,NP),W(NP),V(NP,NP),RV1(NP)
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

