
module helena_boundary
  integer             :: mf, n_bnd
  real, allocatable   :: fr(:), R_bnd(:), Z_bnd(:)
  real                :: Bgeo, Rgeo, Zgeo, amin, eps, ellip, tria_u, tria_l, quad_u, quad_l
  real                :: Reast, Rwest
end

module helena_profiles
  integer             :: n_prof
  real,allocatable    :: psi(:),dp_dpsi(:),fdf_dpsi(:),p_psi(:),f_psi(:),zjz_psi(:),q(:)
  real                :: p_bnd
endmodule

program jorek_to_helena
use constants
use helena_boundary
use helena_profiles

implicit none

real*8  :: R_axis,Z_axis,F0,psi_bnd,psi_axis,psi_xpoint
real*8  :: current, beta_p, beta_t, beta_n
integer :: i

read(5,*) R_axis,Z_axis,F0
read(5,*) psi_bnd,psi_axis,psi_xpoint
read(5,*) n_bnd

allocate(r_bnd(n_bnd),z_bnd(n_bnd))

do i=1,n_bnd
  read(5,*) r_bnd(i),z_bnd(i)
!  write(*,*) r_bnd(i),z_bnd(i)
enddo

read(5,*) amin, Rgeo, Bgeo
read(5,*) current,beta_p,beta_t,beta_n

write(*,*) amin, Rgeo, Bgeo
write(*,*) current,beta_p,beta_t,beta_n

read(5,*) n_prof

n_prof = n_prof+1

allocate(psi(n_prof),dp_dpsi(n_prof),zjz_psi(n_prof),q(n_prof),p_psi(n_prof),f_psi(n_prof),fdf_dpsi(n_prof))

do i=2,n_prof
  read(5,*) psi(i),dp_dpsi(i),zjz_psi(i),q(i)
enddo

write(*,*) ' psi_axis : ',psi_axis
write(*,*) ' psi_bnd  : ',psi_bnd
write(*,*) ' psi_xp   : ',psi_xpoint

write(*,*) ' amin : ',amin
write(*,*) ' Rgeo : ',Rgeo
write(*,*) ' Bgeo : ',Bgeo
write(*,*) ' Current : ',current
write(*,*) ' Beta_p  : ',beta_p
write(*,*) ' Beta_t  : ',beta_t
write(*,*) ' Beta_n  : ',beta_n

psi(1) = psi_axis
dp_dpsi(1) = dp_dpsi(2)
zjz_psi(1) = zjz_psi(2)
q(1)       = q(2)

do i=1,n_prof
  fdf_dpsi(i) = 1. - (psi(i) - psi(1))/(psi(n_prof) - psi(1))
enddo

mf=128

call fshape

Bgeo    = F0 / Rgeo

write(*,*) ' TEST  B     = ',MU_zero, Rgeo**2, dp_dpsi(1), fdf_dpsi(1)

write(*,*)  ' export to HELENA namelist'

open(20,file='helena.nml')
write(20,*) ' &SHAPE IAS = 1, ISHAPE = 2,'
write(20,'(A,i5,A,i5)') '    MFM= ',mf,', MHARM = ',mf/2
do i=1,mf/2
  write(20,'(A,i4,A,e14.6,A,i4,A,e14.6)') '    FM(',2*i-1,')=',fr(2*i-1)/amin,' FM(',2*i,')=',fr(2*i)/amin
enddo
write(20,*) ' &END'
write(20,*) ' &PROFILE '
write(20,*) '    IPAI=11, EPI=1.0, FPI=1.0'
write(20,*) '    IGAM=11'
write(20,*) '    ICUR=11, ECUR=1.0, FCUR=1.0'
write(20,*) '    NPTS = ',n_prof
do i=1,n_prof
  write(20,'(A,i4,A,e12.4,A,i4,A,e12.4,A,i4,A,e12.4,A,i4,A,e12.4)') &
        '    DPR(',i,') = ',dp_dpsi(i),', DF2(',i,') = ',fdf_dpsi(i),', ZJZ(',i,') = ',zjz_psi(i),', QIN(',i,') = ',q(i)
enddo
write(20,*) ' &END'
write(20,*) ' &PHYS '
write(20,'(A,f10.6)') '   EPS   = ',amin/Rgeo
write(20,'(A,f10.6)') '   ALFA  = ',abs(amin**2 * Bgeo / (psi_bnd-psi_axis))
write(20,'(A,f10.6)') '   B     = ',MU_zero * Rgeo**2 * dp_dpsi(1)/fdf_dpsi(1)
write(20,'(A,f10.6)') '   BETAP = ',beta_p
write(20,'(A,e14.6)') '   XIAB  = ',MU_zero * abs(current) / (amin * Bgeo)
write(20,'(A,f10.6)') '   RVAC  = ',Rgeo
write(20,'(A,f10.6)') '   BVAC  = ',Bgeo
write(20,*) ' &END'
write(20,*) ' &NUM'
write(20,*) '    NR    = 51, NP    = 33, NRMAP  = 101,   NPMAP = 129, NCHI = 128'
write(20,*) '    NRCUR = 51, NPCUR = 33, ERRCUR = 1.e-5, NITER = 100, NMESH = 100'
write(20,*) ' &END'
write(20,*) ' &PRI  NPR1=1 &END '
write(20,*) ' &PLOT NPL1=1 &END '
write(20,*) ' &BALL        &END '
close(20)

end

subroutine fshape
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use tr_module
use helena_boundary
use constants

implicit none

real              :: xj, yj, ga
real, allocatable :: THETA(:), GAMMA(:), XV(:),YV(:)
real              :: angle, error, gamm
real, allocatable :: tht_tmp(:),  fr_tmp(:), work(:)
real, allocatable :: tht_sort(:), fr_sort(:), dfr_sort(:)
integer, allocatable :: index_order(:)
real              :: tht, Rbnd_av, ORbnd_av, values(4)
integer           :: m, igrinv, i, j, ishape, ieast(1), iwest(1), n_bnd_short
parameter (error = 1.e-8)

call tr_allocate(fr,1,mf+2,"fr",CAT_GRID)
allocate(theta(mf),gamma(mf),xv(mf),yv(mf))

  write(*,*) ' fshape : (R,Z) set given on ',n_bnd,' points'

  Reast = maxval(R_bnd)
  Rwest = minval(R_bnd)
  ieast = maxloc(R_bnd)
  iwest = minloc(R_bnd)
  Rgeo = (Reast + Rwest) /2.
  Zgeo = (Z_bnd(ieast(1))+Z_bnd(iwest(1)))/2.
  amin = (Reast - Rwest)/2.

  write(*,'(A,3f12.8)') ' Rgeo, Zgeo : ',Rgeo,Zgeo,amin

  call tr_allocate(tht_tmp,1,n_bnd,"tht_tmp",CAT_GRID)
  call tr_allocate(fr_tmp,1,n_bnd,"fr_tmp",CAT_GRID)
  call tr_allocate(work,1,3*n_bnd+6,"work",CAT_GRID)
  call tr_allocate(tht_sort,1,n_bnd+2,"tht_sort",CAT_GRID)
  call tr_allocate(fr_sort,1,n_bnd+2,"fr_sort",CAT_GRID)
  call tr_allocate(dfr_sort,1,n_bnd+2,"dfr_sort",CAT_GRID)
  call tr_allocate(index_order,1,n_bnd+2,"index_order",CAT_GRID)

  do i=1,n_bnd

    tht_tmp(i) = atan2(Z_bnd(i)-Zgeo,R_bnd(i)-Rgeo)

    fr_tmp(i)  = sqrt((R_bnd(i)-Rgeo)**2 + (Z_bnd(i)-Zgeo)**2)

!    write(*,'(i5,2f12.8)') i,tht_tmp(i),fr_tmp(i)

  enddo

  if (abs(tht_tmp(n_bnd) - tht_tmp(1)) .lt. 1.e-6)  n_bnd = n_bnd - 1

  call qsort2(index_order,n_bnd,tht_tmp)

  do i=1,n_bnd
    tht_sort(i) = tht_tmp(index_order(i))
    fr_sort(i)  = fr_tmp(index_order(i))
!    write(*,*) i,tht_sort(i),fr_sort(i)
  enddo

  n_bnd_short = n_bnd
  do i=2,n_bnd

    if ((tht_sort(i) - tht_sort(1)) .gt. 2.*PI) then
      n_bnd_short = i - 1
      exit
    endif

  enddo

  if (abs(tht_sort(n_bnd_short)- tht_sort(1)) .lt. 1.e-6) then
    tht_sort(n_bnd_short) = tht_sort(1) + 2.*PI
    fr_sort(n_bnd_short)  = fr_sort(1)
  else
    tht_sort(n_bnd_short+1) = tht_sort(1) + 2.*PI
    fr_sort(n_bnd_short+1)  = fr_sort(1)
    n_bnd_short = n_bnd_short + 1
  endif

!  write(*,*) ' n_bnd, n_bnd_short : ',n_bnd, n_bnd_short

  call TB15A(n_bnd_short,tht_sort,fr_sort,dfr_sort,work,6)

  do i = 1, mf

   tht = 2.*PI * float(i-1)/float(mf)

   if (tht .lt. tht_sort(1))          tht = tht + 2.*PI
   if (tht .gt. tht_sort(n_bnd_short)) tht = tht - 2.*PI

   call TG02A(0,n_bnd_short,tht_sort,fr_sort,dfr_sort,tht,values)

   fr(i) = values(1)
  enddo

!-------------- FOURIER COEFFICIENTS FRFNUL AND FRF(M) OF FR(J).
call rft2(fr,mf,1)

do m=1,mf
  fr(m) = 2. * fr(m) / float(mf)
enddo
do m=2,mf,2
  fr(m) = - fr(m)
enddo
RETURN
END

SUBROUTINE QSORT2 (ORD,N,A)
!
!==============SORTS THE ARRAY A(I),I=1,2,...,N BY PUTTING THE
!   ASCENDING ORDER VECTOR IN ORD.  THAT IS ASCENDING ORDERED A
!   IS A(ORD(I)),I=1,2,...,N; DESCENDING ORDER A IS A(ORD(N-I+1)),
!   I=1,2,...,N .  THIS SORT RUNS IN TIME PROPORTIONAL TO N LOG N .
!
!
!     ACM QUICKSORT - ALGORITHM #402 - IMPLEMENTED IN FORTRAN BY
!                                 WILLIAM H. VERITY
!                                 COMPUTATION CENTER
!                                 PENNSYLVANIA STATE UNIVERSITY
!                                 UNIVERSITY PARK, PA.  16802
!     With correction to that algorithm.
!
      IMPLICIT INTEGER (A-Z)
!
      DIMENSION ORD(N),POPLST(2,20)
!
!     To sort different input types change the following
!     specification statements; FOR EXAMPLE,  REAL A(N) or
!     CHARACTER *(L) A(N)  for REAL or CHARACTER sorting
!     respectively  similarly for X,XX,Z,ZZ,Y. L is the
!     character length of the elements of A.
!
      REAL A(N)
      REAL X,XX,Z,ZZ,Y
!
      NDEEP=0
      U1=N
      L1=1
      DO 1  I=1,N
    1 ORD(I)=I
    2 IF (U1.GT.L1) GO TO 3
      RETURN
!
    3 L=L1
      U=U1
!
! PART
!
    4 P=L
      Q=U
      X=A(ORD(P))
      Z=A(ORD(Q))
      IF (X.LE.Z) GO TO 5
      Y=X
      X=Z
      Z=Y
      YP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=YP
    5 IF (U-L.LE.1) GO TO 15
      XX=X
      IX=P
      ZZ=Z
      IZ=Q
!
! LEFT
!
    6 P=P+1
      IF (P.GE.Q) GO TO 7
      X=A(ORD(P))
      IF (X.GE.XX) GO TO 8
      GO TO 6
    7 P=Q-1
      GO TO 13
!
! RIGHT
!
    8 Q=Q-1
      IF (Q.LE.P) GO TO 9
      Z=A(ORD(Q))
      IF (Z.LE.ZZ) GO TO 10
      GO TO 8
    9 Q=P
      P=P-1
      Z=X
      X=A(ORD(P))
!
! DIST
!
   10 IF (X.LE.Z) GO TO 11
      Y=X
      X=Z
      Z=Y
      IP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=IP
   11 IF (X.LE.XX) GO TO 12
      XX=X
      IX=P
   12 IF (Z.GE.ZZ) GO TO 6
      ZZ=Z
      IZ=Q
      GO TO 6
!
! OUT
!
   13 CONTINUE
      IF (.NOT.(P.NE.IX.AND.X.NE.XX)) GO TO 14
      IP=ORD(P)
      ORD(P)=ORD(IX)
      ORD(IX)=IP
   14 CONTINUE
      IF (.NOT.(Q.NE.IZ.AND.Z.NE.ZZ)) GO TO 15
      IQ=ORD(Q)
      ORD(Q)=ORD(IZ)
      ORD(IZ)=IQ
   15 CONTINUE
      IF (U-Q.LE.P-L) GO TO 16
      L1=L
      U1=P-1
      L=Q+1
      GO TO 17
   16 U1=U
      L1=Q+1
      U=P-1
   17 CONTINUE
      IF (U1.LE.L1) GO TO 18
!
! START RECURSIVE CALL
!
      NDEEP=NDEEP+1
      POPLST(1,NDEEP)=U
      POPLST(2,NDEEP)=L
      GO TO 3
   18 IF (U.GT.L) GO TO 4
!
! POP BACK UP IN THE RECURSION LIST
!
      IF (NDEEP.EQ.0) GO TO 2
      U=POPLST(1,NDEEP)
      L=POPLST(2,NDEEP)
      NDEEP=NDEEP-1
      GO TO 18
!
! END QSORT
END

SUBROUTINE RFT2(DATA,NR,KR)
!*****************************************************************
! REAL FOURIER TRANSFORM.                                        *
! INPUT:  NR REAL COEFFICIENTS                                   *
!             DATA(1),DATA(1+KR),....,DATA(1+(NR-1)*KR).         *
! OUTPUT: NR/2+1 COMPLEX COEFFICIENTS                            *
!            (DATA(1),      DATA(1+KR))                          *
!            (DATA(1+2*KR), DATA(1+3*KR))                        *
!             .............................                      *
!            (DATA(1+NR*KR),DATA(1+(NR+1)*KR).                   *
! THE CALLING PROGRAM SHOULD HAVE DATA DIMENSIONED WITH AT LEAST *
! (NR+1)*KR+1 ELEMENTS. (I.E., NR+2 IF INCREMENT KR=1).          *
! LASL ROUTINE MAY 75, CALLING FFT2 AND RTRAN2.                  *
!*****************************************************************
implicit none
real    :: DATA(*)
integer :: kr,nr, ktran

CALL FFT2(DATA(1),DATA(KR+1),NR/2,-(KR+KR))
CALL RTRAN2(DATA,NR,KR,1)
RETURN
END

SUBROUTINE RTRAN2(DATA,NR,KR,KTRAN)
!*****************************************************************
! INTERFACE BETWEEN RFT2, RFI2, AND FFT2.                        *
! THE CALLING PROGRAM SHOULD HAVE DATA DIMENSIONED WITH AT LEAST *
! (NR+1)*KR+1 ELEMENTS.                                          *
! LASL ROUTINE MAY 75, CALLED FROM RFT2 AND RFI2.                *
!*****************************************************************
implicit none
real    :: data(*), theta, dc, ds, ws, wc, sumr, difr, sumi, difi
real    :: tr, ti, wca
integer :: nr, kr, ktran, ks, n, nmax, kmax, k, nk

KS=2*KR
N=NR/2
NMAX=N*KS+2
KMAX=NMAX/2
THETA=1.5707963267949/FLOAT(N)
DC=2.*SIN(THETA)**2
DS=SIN(2.*THETA)
WS=0.

IF (KTRAN .LE. 0) THEN
   WC=-1.0
   DS=-DS
ELSE
   WC=1.0
   DATA(NMAX-1)=DATA(1)
   DATA(NMAX-1+KR)=DATA(KR+1)
ENDIF
DO K=1,KMAX,KS
   NK=NMAX-K
   SUMR=.5*(DATA(K)+DATA(NK))
   DIFR=.5*(DATA(K)-DATA(NK))
   SUMI=.5*(DATA(K+KR)+DATA(NK+KR))
   DIFI=.5*(DATA(K+KR)-DATA(NK+KR))
   TR=WC*SUMI-WS*DIFR
   TI=WS*SUMI+WC*DIFR
   DATA(K)=SUMR+TR
   DATA(K+KR)=DIFI-TI
   DATA(NK)=SUMR-TR
   DATA(NK+KR)=-DIFI-TI
   WCA=WC-DC*WC-DS*WS
   WS=WS+DS*WC-DC*WS
   WC=WCA
enddo
return
end

SUBROUTINE FFT2 (DATAR,DATAI,N,INC)
!*****************************************************************
! FFT2 FORTRAN VERSION CLAIR NIELSON MAY 75.                     *
!*****************************************************************
real    :: DATAR(*), DATAI(*)
integer :: n, ninc

KTRAN=ISIGN(-1,INC)
KS=IABS(INC)
IP0=KS
IP3=IP0*N
IREV=1

      DO I=1,IP3,IP0
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
         IREV=IREV+IBIT
      enddo
      IP1=IP0
      THETA=REAL(KTRAN)*3.1415926535898
   30 IF(IP1.GE.IP3) return
      IP2=IP1+IP1
      SINTH=SIN(.5*THETA)
      WSTPR=-2.*SINTH*SINTH
      WSTPI=SIN(THETA)
      WR=1.
      WI=0.
      DO I1=1,IP1,IP0
         DO I3=I1,IP3,IP2
            J0=I3
            J1=J0+IP1
            TEMPR=WR*DATAR(J1)-WI*DATAI(J1)
            TEMPI=WR*DATAI(J1)+WI*DATAR(J1)
            DATAR(J1)=DATAR(J0)-TEMPR
            DATAI(J1)=DATAI(J0)-TEMPI
            DATAR(J0)=DATAR(J0)+TEMPR
            DATAI(J0)=DATAI(J0)+TEMPI
         enddo
         TEMPR=WR
         WR=WR*WSTPR-WI*WSTPI+WR
         WI=WI*WSTPR+TEMPR*WSTPI+WI
      enddo
      IP1=IP2
      THETA=.5*THETA
      GOTO 30
RETURN
END


SUBROUTINE FSUM2(F,T,FFNUL,FFCOS,FFSIN,MHARM)
!-----------------------------------------------------------------------
! FOURIER SYNTHESIS OF GENERAL  FUNCTION F(T) AT SINGLE POINT T.
!-----------------------------------------------------------------------
implicit none
integer :: mharm, m
real    :: ffnul, ffcos(*), ffsin(*), f, t, s, c, co, ca, si, sum

CO=COS(T)
SI=SIN(T)
C=1.
S=0.
SUM=.5*FFNUL
do m=1,mharm
  CA=C*CO-S*SI
  S=S*CO+C*SI
  C=CA
  SUM=SUM+FFCOS(M)*C + FFSIN(M)*S
enddo
F=SUM
RETURN
END


SUBROUTINE TB15A(N,X,F,D,W,LP)
!------------------------------------------------------------------
! HSL routine for cubic spline with periodic boundary conditions
! first point must be the same as last : f(1)=f(n)
!    N : number of points
!    X : coordinate (input)
!    F : the function values to be splined (input)
!    D : the derivatives at the points (output)
!    W : workspace (dimension 3N)
!   LP : unit number for output
!------------------------------------------------------------------
REAL ZERO,ONE,TWO,THREE
PARAMETER (ZERO=0.0E0,ONE=1.0E0,TWO=2.0E0,THREE=3.0E0)
INTEGER LP,N
REAL D(N),F(N),W(*),X(N)
REAL A3N1,F1,F2,H1,H2,P
INTEGER I,J,K,N2

WRITE(*,*) F(1),F(N)
IF (N.LT.4) THEN
  WRITE (LP,'(A39)')  'RETURN FROM TB15AD BECAUSE N TOO SMALL'
  W(1) = ONE
  RETURN
END IF
DO I = 2,N
  IF (X(I).LE.X(I-1)) THEN
    WRITE (LP,'(A29,I3,A13)') ' RETURN FROM TB15AD BECAUSE  ',I,' OUT OF ORDER'
    W(1) = TWO
    RETURN
  END IF
ENDDO
IF (F(1).NE.F(N)) THEN
  WRITE (LP,'(A40)')  'RETURN FROM TB15AD BECAUSE F(1).NE.F(N)'
  W(1) = THREE
  RETURN
END IF
DO I = 2,N
  H1 = ONE/ (X(I)-X(I-1))
  F1 = F(I-1)
  IF (I.EQ.N) THEN
    H2 = ONE/ (X(2)-X(1))
    F2 = F(2)
  ELSE
    H2 = ONE/ (X(I+1)-X(I))
    F2 = F(I+1)
  END IF
  W(3*I-2) = H1
  W(3*I-1) = TWO* (H1+H2)
  W(3*I) = H2
  D(I) = 3.0* (F2*H2*H2+F(I)* (H1*H1-H2*H2)-F1*H1*H1)
ENDDO
N2 = N - 2
K = 5
A3N1 = W(3*N-1)
DO I = 2,N2
  P = W(K+2)/W(K)
  W(K+3) = W(K+3) - P*W(K+1)
  D(I+1) = D(I+1) - P*D(I)
  W(K+2) = -P*W(K-1)
  P = W(K-1)/W(K)
  A3N1 = -P*W(K-1) + A3N1
  D(N) = D(N) - P*D(I)
  K = K + 3
ENDDO
P = (W(K+2)+W(K-1))/W(K)
A3N1 = A3N1 - P* (W(K+1)+W(K-1))
D(N) = (D(N)-P*D(N-1))/A3N1
DO I = 3,N
  J = N + 2 - I
  D(J) = (D(J)-W(3*J)*D(J+1)-W(3*J-2)*D(N))/W(3*J-1)
ENDDO
D(1) = D(N)
W(1) = ZERO
RETURN
END



SUBROUTINE TG02A(IX,N,U,S,D,X,V)
!------------------------------------------------------------------
! HSL subroutine to calculate splined values
!    N  : number of points
!    IX : negative 0 -> no initial guess for where xi is
!        positive -> gues for index close to value X
!   U   : the coordinates of the spline points
!   S   : the function values of the spline points
!   D   : the derivatives on the spline points
!   X   : the coordinate where the output is wanted
!   V(1-4) : value and derivatives of the spline interpolation
!------------------------------------------------------------------
REAL X
INTEGER IX,N
REAL D(*),S(*),U(*),V(*)
REAL A,B,C,C3,D0,D1,EPS,GAMA,H,HR,HRR,PHI,S0,S1,T,THETA
INTEGER I,IFLG,J

EPS = 1.E-33
K = 0
IFLG = 0
IF (X.LT.U(1)) GO TO 990
IF (X.GT.U(N)) GO TO 991
IF (IX.LT.0 .OR. IFLG.EQ.0) GO TO 12
IF (X.GT.U(J+1)) GO TO 1
IF (X.GE.U(J)) GO TO 18
GO TO 2

    1 J = J + 1
   11 IF (X.GT.U(J+1)) GO TO 1
      GO TO 7
   12 J = ABS(X-U(1))/ (U(N)-U(1))* (N-1) + 1
      J = MIN(J,N-1)
      IFLG = 1
      IF (X.GE.U(J)) GO TO 11
    2 J = J - 1
      IF (X.LT.U(J)) GO TO 2
    7 K = J
      H = U(J+1) - U(J)
      HR = 1./H
      HRR = (HR+HR)*HR
      S0 = S(J)
      S1 = S(J+1)
      D0 = D(J)
      D1 = D(J+1)
      A = S1 - S0
      B = A - H*D1
      A = A - H*D0
      C = A + B
      C3 = C*3.
   18 THETA = (X-U(J))*HR
      PHI = 1. - THETA
      T = THETA*PHI
      GAMA = THETA*B - PHI*A
      V(1) = THETA*S1 + PHI*S0 + T*GAMA
      V(2) = THETA*D1 + PHI*D0 + T*C3*HR
      V(3) = (C* (PHI-THETA)-GAMA)*HRR
      V(4) = -C3*HRR*HR
      RETURN
  990 IF (X.LE.U(1)-EPS*MAX(ABS(U(1)),ABS(U(N)))) GO TO 99
      J = 1
      GO TO 7
  991 IF (X.GE.U(N)+EPS*MAX(ABS(U(1)),ABS(U(N)))) GO TO 995
      J = N - 1
      GO TO 7
  995 K = N
   99 IFLG = 0
      DO I = 1,4
        V(I) = 0.
      ENDDO
RETURN
END
