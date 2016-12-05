      SUBROUTINE RFT2(DATA,NR,KR)
      !*****************************************************************
      !* REAL FOURIER TRANSFORM.                                       *
      !* INPUT:  NR REAL COEFFICIENTS                                  *
      !*            DATA(1),DATA(1+KR),....,DATA(1+(NR-1)*KR).         *
      !* OUTPUT: NR/2+1 COMPLEX COEFFICIENTS                           *
      !*           (DATA(1),      DATA(1+KR))                          *
      !*           (DATA(1+2*KR), DATA(1+3*KR))                        *
      !*            .............................                      *
      !*           (DATA(1+NR*KR),DATA(1+(NR+1)*KR).                   *
      !* THE CALLING PROGRAM SHOULD HAVE DATA DIMENSIONED WITH AT LEAST*
      !* (NR+1)*KR+1 ELEMENTS. (I.E., NR+2 IF INCREMENT KR=1).         *
      !* LASL ROUTINE MAY 75, CALLING FFT2 AND RTRAN2.                 *
      !*****************************************************************

      real*8  :: DATA(*)
      integer :: nr, kr

      CALL FFT2(DATA(1),DATA(KR+1),NR/2,-(KR+KR))
      CALL RTRAN2(DATA,NR,KR,1)

      RETURN
      END

      SUBROUTINE RTRAN2(DATA,NR,KR,KTRAN)
      !*****************************************************************
      !* INTERFACE BETWEEN RFT2, RFI2, AND FFT2.                       *
      !* THE CALLING PROGRAM SHOULD HAVE DATA DIMENSIONED WITH AT LEAST*
      !* (NR+1)*KR+1 ELEMENTS.                                         *
      !* LASL ROUTINE MAY 75, CALLED FROM RFT2 AND RFI2.               *
      !*****************************************************************

      real*8 DATA(*)

      PI = 2.d0*asin(1.d0)

      KS=2*KR
      N=NR/2
      NMAX=N*KS+2
      KMAX=NMAX/2
      THETA=PI/FLOAT(2*N)
      DC=2.D0*SIN(THETA)**2
      DS=SIN(2.D0*THETA)
      WS=0.D0
      IF(KTRAN.LE.0) THEN
         WC=-1.0D0
         DS=-DS
      ELSE
         WC=1.D0
         DATA(NMAX-1)=DATA(1)
         DATA(NMAX-1+KR)=DATA(KR+1)
      ENDIF
      DO K=1,KMAX,KS
         NK=NMAX-K
         SUMR=.5D0*(DATA(K)+DATA(NK))
         DIFR=.5D0*(DATA(K)-DATA(NK))
         SUMI=.5D0*(DATA(K+KR)+DATA(NK+KR))
         DIFI=.5D0*(DATA(K+KR)-DATA(NK+KR))
         TR=WC*SUMI-WS*DIFR
         TI=WS*SUMI+WC*DIFR
         DATA(K)=SUMR+TR
         DATA(K+KR)=DIFI-TI
         DATA(NK)=SUMR-TR
         DATA(NK+KR)=-DIFI-TI
         WCA=WC-DC*WC-DS*WS
         WS=WS+DS*WC-DC*WS
         WC=WCA
      ENDDO
      RETURN
      END

      SUBROUTINE FFT2 (DATAR,DATAI,N,INC)
      !*****************************************************************
      !* FFT2 FORTRAN VERSION CLAIR NIELSON MAY 75.                    *
      !*****************************************************************

      DIMENSION DATAR(*), DATAI(*)

      PI = 2.d0*asin(1.d0)

      KTRAN=ISIGN(-1,INC)
      KS=IABS(INC)
      IP0=KS
      IP3=IP0*N
      IREV=1
      DO 20 I=1,IP3,IP0
         IF(I.LT.IREV) THEN
            TEMPR=DATAR(I)
            TEMPI=DATAI(I)
            DATAR(I)=DATAR(IREV)
            DATAI(I)=DATAI(IREV)
            DATAR(IREV)=TEMPR
            DATAI(IREV)=TEMPI
         ENDIF
         IBIT=IP3/2
   10    IF(IREV.GT.IBIT) THEN
            IREV=IREV-IBIT
            IBIT=IBIT/2
            IF(IBIT.GE.IP0) GOTO 10
         ENDIF
   20    IREV=IREV+IBIT
      IP1=IP0
      THETA=REAL(KTRAN)*PI
   30 IF(IP1.GE.IP3) GOTO 60
      IP2=IP1+IP1
      SINTH=SIN(.5D0*THETA)
      WSTPR=-2.D0*SINTH*SINTH
      WSTPI=SIN(THETA)
      WR=1.D0
      WI=0.D0
      DO 50 I1=1,IP1,IP0
         DO 40 I3=I1,IP3,IP2
            J0=I3
            J1=J0+IP1
            TEMPR=WR*DATAR(J1)-WI*DATAI(J1)
            TEMPI=WR*DATAI(J1)+WI*DATAR(J1)
            DATAR(J1)=DATAR(J0)-TEMPR
            DATAI(J1)=DATAI(J0)-TEMPI
            DATAR(J0)=DATAR(J0)+TEMPR
   40       DATAI(J0)=DATAI(J0)+TEMPI
         TEMPR=WR
         WR=WR*WSTPR-WI*WSTPI+WR
   50    WI=WI*WSTPR+TEMPR*WSTPI+WI
      IP1=IP2
      THETA=.5D0*THETA
      GOTO 30
   60 RETURN

      END
