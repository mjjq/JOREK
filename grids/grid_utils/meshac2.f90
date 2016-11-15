SUBROUTINE MESHAC2(NR,SG,XR1,XR2,SIG1,SIG2,BGF,FACT)
!-----------------------------------------------------------------------
! subroutine to construct non-equidistant radial mesh
! the starting value is given by SG(1)
!-----------------------------------------------------------------------
use fgauss
implicit none
integer :: NMAX
parameter (NMAX = 10001)

real*8  :: SG(*), XR1, XR2, SIG1, SIG2, BGF, FACT, SUM, FINT2, FINT1, DSI, DFG, FI
real*8  :: S1(NMAX),FSUM(NMAX)
integer :: I, J, NR

!write(*,*) '*******************************'
!write(*,*) '* MESHAC                      *'
!write(*,'(4f16.6)') xr1,sig1,xr2,sig2
!write(*,*) '*******************************'


if (NR .le. 0) return    ! do nothing for zero points

!--------------------------------------- set parameters of gaussians
!BGF  = 0.6d0
!FACT = 1.d0

if ((SG(1).LT. 0.d0) .OR. (SG(1) .GE. 1.d0)) SG(1) = 0.d0

!--------------------------------------- integrate gaussian
DSI      = (1.d0 - SG(1)) / FLOAT(NMAX-1)
S1(1)    = SG(1)
FSUM(1)  = 0.d0
SUM      = 0.d0
FINT2 = FGAUS(S1(1),BGF,XR1,XR2,SIG1,SIG2,FACT,DFG)
DO I=2,NMAX
  S1(I) = SG(1) + FLOAT(I-1) * DSI
  FINT1 = FINT2
  FINT2 = FGAUS(S1(I),BGF,XR1,XR2,SIG1,SIG2,FACT,DFG)
  SUM = SUM + (FINT1+FINT2)/2.d0 * DSI
  FSUM(I) = SUM
ENDDO

DO I=1,NMAX-1
  FSUM(I) = FSUM(I)/FSUM(NMAX)
ENDDO
FSUM(NMAX) = 1.d0

J = 2
DO I=2,NR-1
  FI = FLOAT(I-1)/FLOAT(NR-1)
  do while ((FI .gt. FSUM(J)) .and. (J .lt. NMAX))
    J = J + 1
  enddo
  SG(I)   = S1(J-1) + (FI-FSUM(J-1))/(FSUM(J)-FSUM(J-1))*(S1(J)-S1(J-1))
ENDDO
SG(NR)   = 1.d0
SG(1)    = 0.d0

RETURN
END

SUBROUTINE MESHAC3(NR,SG,XR1,XR2,XR3,SIG1,SIG2,SIG3,BGF,FACT)
!-----------------------------------------------------------------------
! subroutine to construct non-equidistant radial mesh
! the starting value is given by SG(1)
!-----------------------------------------------------------------------
use fgauss
implicit none
integer :: NMAX
parameter (NMAX = 10001)

real*8  :: SG(*), XR1, XR2, XR3, SIG1, SIG2, SIG3
real*8  :: BGF, FACT, SUM, FINT2, FINT1, DSI, DFG, FI
real*8  :: S1(NMAX),FSUM(NMAX)
integer :: I, J, NR

write(*,*) '*******************************'
write(*,*) '* MESHAC3                      *'
write(*,'(6f16.6)') xr1,sig1,xr2,sig2,xr3,sig3
write(*,*) '*******************************'


if (NR .le. 0) return    ! do nothing for zero points

!--------------------------------------- set parameters of gaussians
!BGF  = 0.6d0
!FACT = 1.d0

if ((SG(1).LT. 0.d0) .OR. (SG(1) .GE. 1.d0)) SG(1) = 0.d0

!--------------------------------------- integrate gaussian
DSI      = (1.d0 - SG(1)) / FLOAT(NMAX-1)
S1(1)    = SG(1)
FSUM(1)  = 0.d0
SUM      = 0.d0
FINT2 = FGAUS3(S1(1),BGF,XR1,XR2,XR3,SIG1,SIG2,SIG3,FACT,DFG)
DO I=2,NMAX
  S1(I) = SG(1) + FLOAT(I-1) * DSI
  FINT1 = FINT2
  FINT2 = FGAUS3(S1(I),BGF,XR1,XR2,XR3,SIG1,SIG2,SIG3,FACT,DFG)
  SUM = SUM + (FINT1+FINT2)/2.d0 * DSI
  FSUM(I) = SUM
ENDDO

DO I=1,NMAX-1
  FSUM(I) = FSUM(I)/FSUM(NMAX)
ENDDO
FSUM(NMAX) = 1.d0

J = 2
DO I=2,NR-1
  FI = FLOAT(I-1)/FLOAT(NR-1)
  do while ((FI .gt. FSUM(J)) .and. (J .lt. NMAX))
    J = J + 1
  enddo
  SG(I)   = S1(J-1) + (FI-FSUM(J-1))/(FSUM(J)-FSUM(J-1))*(S1(J)-S1(J-1))
ENDDO
SG(NR)   = 1.d0
SG(1)    = 0.d0

RETURN
END
